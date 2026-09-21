r"""The melt callbacks float ice with the density the forward floats it with.

Until September 2026 they tested flotation with the fresh-water density of the
SMB unit conversion, which put the flotation surface 2.4 percent of the water
depth too low and read every floating cell thicker than 78 percent of its
flotation thickness as grounded, withholding its melt (icepack/ismip7#66).
"""
import numpy as np
import pytest

from icepack2_tools.forcing import height_above_flotation, is_floating

RHO_I, RHO_SW = 917.0, 1024.0


def _surface(b, h):
    r"""The forward's surface: grounded on the bed, or hydrostatic."""
    return np.maximum(b + h, (1.0 - RHO_I / RHO_SW) * h)


def test_ice_at_its_flotation_thickness_floats_and_thicker_ice_grounds():
    b = np.array([-500.0, -500.0, -500.0])
    h_float = -b * RHO_SW / RHO_I                 # exactly at flotation
    h = np.array([h_float[0] * 0.5, h_float[1] * 0.99, h_float[2] * 1.01])
    afloat = is_floating(_surface(b, h), b)
    assert afloat.tolist() == [True, True, False]


def test_the_grounding_zone_band_melts():
    r"""A cell at 0.9 times its flotation thickness floats, and the forward
    floats it; the fresh-water test read it as grounded, so it got no melt."""
    b = np.array([-570.0])
    h = 0.9 * (-b) * RHO_SW / RHO_I
    s = _surface(b, h)
    assert is_floating(s, b).all()
    assert height_above_flotation(s, b, rho_water=1000.0)[0] > 0.0


def test_ice_free_ocean_reads_as_floating():
    r"""A cell holding no ice over deep water has haf < 0; the callers that
    need ice present test h > 0 themselves."""
    assert is_floating(_surface(np.array([-500.0]), np.array([0.0])), np.array([-500.0])).all()


def test_grounded_ice_above_sea_level_does_not_float():
    b = np.array([100.0])
    assert not is_floating(_surface(b, np.array([50.0])), b).any()
