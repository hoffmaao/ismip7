r"""The front-flux probe's own arithmetic.

``probe_front_flux.py`` answers whether a calving law acts on ice that
matters, so its numbers have to be the ones the level set would remove: the
flux-weighted front thickness, the outward normal speed, and
``sum max(u.n, 0) h L`` over the front cells. Checked on an imposed state
where every number follows by hand.
"""
import importlib.util
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_probe():
    r"""Loaded at import time, before any mesh exists: the module's own import
    chain reaches Irksome, which refuses to be imported after a UFL
    multifunction has run."""
    path = os.path.join(REPO, "antarctica", "scripts", "probe_front_flux.py")
    spec = importlib.util.spec_from_file_location("probe_front_flux", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


probe = _load_probe()

import numpy as np  # noqa: E402
import pytest  # noqa: E402
from firedrake import (Constant, Function, FunctionSpace,  # noqa: E402
                       SpatialCoordinate, UnitSquareMesh, VectorFunctionSpace,
                       as_vector)


@pytest.fixture(scope="module")
def state():
    r"""Ice on the left half of a unit domain, 300 m thick, flowing out at
    1 m/yr, so the front is the column at x = 0.5."""
    mesh = UnitSquareMesh(8, 8)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)
    h = Function(Q0).interpolate(300.0 * (x < 0.5))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((Constant(1.0), Constant(0.0))))
    return mesh, h, u


def test_the_front_is_found_and_its_thickness_is_the_ice_thickness(state):
    m = probe
    mesh, h, u = state
    cells, length, h_front, u_front, flux = m.front_flux(mesh, h, u, 1.0)
    assert cells > 0
    assert length > 0.0
    assert h_front == pytest.approx(300.0)
    assert flux > 0.0


def test_the_flux_is_thickness_times_speed_times_front_length(state):
    r"""With a uniform outward speed and thickness the flux is the product,
    which is what makes it comparable with an observed calving flux."""
    m = probe
    mesh, h, u = state
    _, length, h_front, u_front, flux = m.front_flux(mesh, h, u, 1.0)
    # u . n at a front whose normal is +x and whose speed is 1 m/yr
    assert u_front == pytest.approx(1.0, rel=1e-6)
    assert flux == pytest.approx(h_front * u_front * length * m.RHO_GT, rel=1e-6)


def test_a_threshold_above_the_ice_leaves_no_front(state):
    m = probe
    mesh, h, u = state
    cells, length, _, _, flux = m.front_flux(mesh, h, u, 400.0)
    assert cells == 0
    assert flux == 0.0


def test_inflow_does_not_count_as_calving(state):
    r"""Only ice leaving counts: a front the ice flows INTO removes nothing,
    and a probe that summed the signed flux would report a negative calving
    flux and hide the real one."""
    m = probe
    mesh, h, _ = state
    u_in = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((Constant(-1.0), Constant(0.0))))
    _, _, _, u_front, flux = m.front_flux(mesh, h, u_in, 1.0)
    assert u_front < 0.0
    assert flux == 0.0
