r"""The ISMIP7 SMB anomaly enters as the mass flux over the ice density.

``acabf`` and ``acabf-anomaly`` are in kg m-2 s-1. The thickness tendency the
transport needs is that flux over the ice density, 917 kg m-3, which is also
how ``write_ismip7_output`` converts the model's ``acabf`` back for submission.
Until September 2026 the reader multiplied by a further water-to-ice density
ratio, so every ESM anomaly entered 9 percent too large (issue #29).
"""
import pytest

from icepack2_tools.forcing import _RHO_ICE, _SEC_PER_YEAR, smb_kgm2s_to_myr


def test_one_metre_of_ice_per_year_round_trips():
    flux = _RHO_ICE / _SEC_PER_YEAR              # the flux that is 1 m/yr of ice
    assert smb_kgm2s_to_myr(flux) == pytest.approx(1.0, rel=1e-12)


def test_the_conversion_is_the_flux_over_the_ice_density():
    assert smb_kgm2s_to_myr(1.0) == pytest.approx(_SEC_PER_YEAR / _RHO_ICE, rel=1e-12)
    inflated = _SEC_PER_YEAR / _RHO_ICE * (1000.0 / 917.0)   # the pre-fix reading
    assert smb_kgm2s_to_myr(1.0) != pytest.approx(inflated, rel=1e-3)


def test_the_output_writer_inverts_it_up_to_the_calendar():
    out = pytest.importorskip("icepack2_tools.ismip7_output")
    back = out.RHO_I / out.SECONDS_PER_YEAR      # m ice/yr -> kg m-2 s-1
    # 31556926 s in, 31557600 s out: the 2e-5 the submission README states.
    assert smb_kgm2s_to_myr(0.01) * back == pytest.approx(0.01, rel=1e-4)
