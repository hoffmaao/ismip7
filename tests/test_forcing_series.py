r"""The atmosphere reader finds the renamed MRI product and persists the last
forcing year past the end of the series.

Discussion #37 (Aug 2026) renamed MRI-ESM2-0's downscaled atmosphere from
``SDBN1-<res>`` to ``GEMB-SDBN1-<res>`` (data unchanged); discussion #8:
CESM2-WACCM stops at 2299 and the empty 2300 files were removed, while a
2015-2300 run needs the 2300 forcing year, so exactly one year past the end
of the series is bridged and anything further is a short tree and an error.
All of it is exercised on a three-by-three synthetic tree written with
xarray.
"""
import os

import numpy as np
import pytest

xr = pytest.importorskip("xarray")

from icepack2_tools.forcing import (ISMIP7Atmosphere, atmosphere_product,  # noqa: E402
                                    forcing_year)


def _write_year(vdir, esm, scenario, product, version, year, value):
    times = np.array([np.datetime64(f"{year}-{m:02d}-15") for m in range(1, 13)])
    ds = xr.Dataset(
        {"acabf": (("time", "y", "x"), np.full((12, 3, 3), value, dtype="float32"))},
        coords={"time": times, "y": np.array([0.0, 8000.0, 16000.0]), "x": np.array([0.0, 8000.0, 16000.0])},
    )
    ds.to_netcdf(os.path.join(vdir, f"acabf_AIS_{esm}_{scenario}_{product}_{version}_{year}.nc"))


@pytest.fixture
def mri_tree(tmp_path):
    esm, scenario, product, version = "MRI-ESM2-0", "ssp585", "GEMB-SDBN1-8000m", "v1"
    vdir = tmp_path / esm / scenario / product / "acabf" / version
    vdir.mkdir(parents=True)
    _write_year(str(vdir), esm, scenario, product, version, 2298, 1.0)
    _write_year(str(vdir), esm, scenario, product, version, 2299, 2.0)
    return tmp_path


def test_renamed_product_is_found(mri_tree):
    assert atmosphere_product(str(mri_tree), "MRI-ESM2-0", "ssp585") == "GEMB-SDBN1-8000m"
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    field = atm._load_year("acabf", 2299)
    assert field is not None and np.allclose(np.asarray(field), 2.0)


def test_last_year_persists_past_the_series(mri_tree, capsys):
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    f2300 = atm._load_year("acabf", 2300)
    f2299 = atm._load_year("acabf", 2299)
    assert f2300 is not None and np.allclose(np.asarray(f2300), np.asarray(f2299))
    out = capsys.readouterr().out
    assert out.count("persisting 2299") == 1, out           # reported once per variable


def test_further_past_the_series_is_an_error(mri_tree):
    r"""The bridge is exactly one year: a tree that stopped short must fail
    rather than repeat its last year for decades."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    with pytest.raises(FileNotFoundError, match=r"acabf.*no year 2305"):
        atm._load_year("acabf", 2305)


def test_a_gap_inside_the_series_is_an_error(mri_tree):
    r"""A hole in a series that exists must never reach get_field, which
    would hand the run a field of zeros for that year."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    with pytest.raises(FileNotFoundError, match=r"acabf.*no year 2200"):
        atm._load_year("acabf", 2200)
    with pytest.raises(FileNotFoundError):
        atm.get_smb(2200, np.zeros(3), np.zeros(3), anomaly=False)


def test_a_year_before_the_series_is_reported_as_such(mri_tree):
    r"""A projection asking its scenario tree for a year the scenario does not
    start until is a timeline error, not a short download, and the message
    has to say which so the operator fixes the right thing."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    with pytest.raises(FileNotFoundError, match=r"no year 2014.*precedes the series"):
        atm._load_year("acabf", 2014)


def test_an_absent_variable_stays_optional(mri_tree):
    r"""A variable with no files at all is an optional product (dacabfdz,
    ts-anomaly), not a gap, so it still reads as absent."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    assert atm._load_year("dacabfdz", 2299) is None


# ---- which year a step's forcing comes from ------------------------------
#
# run_simulation hands the callback the END of the step, so the step from
# 2300.9 to 2301.0 lies in 2300. Rounding instead pushed every step past the
# half-year mark into the next year, which made a 2015-2300 projection end by
# requesting 2301: two past CESM2-WACCM's 2299, so the one-year end-of-series
# bridge could not cover it and the run died on its last step.


@pytest.mark.parametrize("t_yr,expected", [
    (2015.1, 2015),          # the first step of a run starting at 2015.0
    (2015.5, 2015),
    (2015.9, 2015),
    (2016.0, 2015),          # a step ENDING on 1 January belongs to the year before
    (2016.1, 2016),
    (2300.6, 2300),          # used to round up to 2301
    (2301.0, 2300),          # the last step of a 2015-2300 projection
])
def test_the_forcing_year_is_the_year_the_step_lies_in(t_yr, expected):
    assert forcing_year(t_yr) == expected


def test_every_step_of_a_year_draws_the_same_forcing_year():
    r"""run_simulation steps t_yr = t_start + k*dt for k >= 1, so a year's ten
    dt=0.1 steps end at Y.1 .. Y+1.0 and all ten belong to year Y. Rounding
    split them five and five across two forcing years."""
    steps = [round(2015.0 + k * 0.1, 10) for k in range(1, 11)]
    assert steps[-1] == 2016.0
    assert {forcing_year(t) for t in steps} == {2015}


def test_the_last_step_of_a_cesm_projection_reaches_the_bridge(mri_tree):
    r"""A projection covering 2015-2300 ends at t=2301.0. Against a series
    that stops at 2299 that has to land on year 2300, which the one-year
    persistence covers; 2301 would be two past and would raise."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    field = atm.get_field("acabf", 2301.0, np.zeros(3), np.zeros(3))
    assert np.allclose(field, 2.0)         # the 2299 file, persisted for 2300
