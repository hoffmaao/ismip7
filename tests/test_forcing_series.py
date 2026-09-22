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
                                    forcing_year, make_forcing_callback)


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
    # 2296, then a HOLE at 2297, then 2298-2299: the span has an interior gap
    # so the "missing from the series" branch can be probed for what it is
    _write_year(str(vdir), esm, scenario, product, version, 2296, 0.5)
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
    would hand the run a field of zeros for that year. 2297 is INSIDE the
    fixture's 2296-2299 span, so this is the missing-year branch and not the
    precedes-the-series one."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    with pytest.raises(FileNotFoundError, match=r"no year 2297.*missing from the series"):
        atm._load_year("acabf", 2297)
    with pytest.raises(FileNotFoundError, match="missing from the series"):
        atm.get_smb(2297, np.zeros(3), np.zeros(3), anomaly=False)


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


def test_get_smb_divides_the_flux_by_the_ice_density(mri_tree):
    r"""``acabf`` is a mass flux, so the field the transport receives is the
    flux over the ice density and nothing else (issue #29)."""
    from icepack2_tools.forcing import _RHO_ICE, _SEC_PER_YEAR
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    smb = atm.get_smb(2299, np.zeros(3), np.zeros(3), anomaly=False)
    assert np.allclose(smb, 2.0 * _SEC_PER_YEAR / _RHO_ICE, rtol=1e-12)


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


def test_the_reader_takes_a_plain_calendar_year(mri_tree):
    r"""The readers are also called with years straight out of
    ``available_years`` (the climatology pools in ``smb_scheme`` and
    ``compute_climatology``). Shifting inside the reader re-referenced every
    pooled year by one and made a scenario's first year read as preceding its
    own series, which those callers swallow into a silent SMB fallback."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    assert np.allclose(atm.get_field("acabf", 2296, np.zeros(3), np.zeros(3)), 0.5)
    assert np.allclose(atm.get_field("acabf", 2298, np.zeros(3), np.zeros(3)), 1.0)
    assert np.allclose(atm.get_field("acabf", 2299, np.zeros(3), np.zeros(3)), 2.0)
    # the first year of the series resolves to itself, not to the year before
    assert np.allclose(atm.get_smb(2296, np.zeros(3), np.zeros(3), anomaly=False),
                       atm.get_smb(2296.0, np.zeros(3), np.zeros(3), anomaly=False))


def test_the_last_step_of_a_cesm_projection_reaches_the_bridge(mri_tree):
    r"""A projection covering 2015-2300 ends at t=2301.0. Composed the way the
    forcing callback composes them, that model time lands on year 2300, which
    the one-year persistence covers; asking the reader for 2301 outright would
    be two past the series and would raise."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    field = atm.get_field("acabf", forcing_year(2301.0), np.zeros(3), np.zeros(3))
    assert np.allclose(field, 2.0)         # the 2299 file, persisted for 2300
    with pytest.raises(FileNotFoundError):
        atm.get_field("acabf", 2301, np.zeros(3), np.zeros(3))


class _RecordingAtmosphere:
    r"""Stands in for ISMIP7Atmosphere and records the year it was asked for."""

    def __init__(self):
        self.asked = []

    def get_smb(self, year, mesh_x, mesh_y, anomaly=True):
        self.asked.append(year)
        return np.zeros(len(mesh_x))


class _Accum:
    class _Dat:
        def __init__(self, n):
            self.data = np.zeros(n)

    def __init__(self, n):
        self.dat = self._Dat(n)


def test_the_callback_is_the_only_place_model_time_becomes_a_year():
    r"""The conversion belongs at the boundary that knows it is holding model
    time. Driving the real callback over one year of dt=0.1 steps, every step
    asks for that one calendar year, and the last step of a 2015-2300 run asks
    for 2300 rather than 2301."""
    atm = _RecordingAtmosphere()
    callback = make_forcing_callback(atm=atm)
    ctx = {"geom_xy": (np.zeros(3), np.zeros(3)), "accum": _Accum(3)}

    for k in range(1, 11):
        callback(ctx, round(2015.0 + k * 0.1, 10))
    assert atm.asked == [2015] * 10
    assert all(isinstance(y, int) for y in atm.asked)

    atm.asked.clear()
    callback(ctx, 2301.0)
    assert atm.asked == [2300]


# --- what a run opened, for the committed report ---------------------------

def test_provenance_names_what_the_reader_would_open(mri_tree):
    from icepack2_tools.forcing import FORCING_PROVENANCE_MARKER, describe_forcing_provenance
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585")
    # the reader asks for v2 first; only v1 is on disk, so v1 is what it opens
    (row,) = atm.provenance(("acabf", "acabf-anomaly"))
    assert (row["variable"], row["product"], row["version"], row["newer"]) == ("acabf", "GEMB-SDBN1-8000m", "v1", [])
    (line,) = describe_forcing_provenance(atm, None, variables={"atmosphere": ("acabf",)})
    assert line == f"{FORCING_PROVENANCE_MARKER} atmosphere acabf MRI-ESM2-0 ssp585 GEMB-SDBN1-8000m v1"


def test_a_newer_version_beside_the_pin_is_said_out_loud(mri_tree):
    from icepack2_tools.forcing import describe_forcing_provenance
    parent = mri_tree / "MRI-ESM2-0" / "ssp585" / "GEMB-SDBN1-8000m" / "acabf"
    (parent / "v2").mkdir()
    (parent / "v3").mkdir()
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585")
    # the pin holds, so a campaign does not change forcing under itself...
    (row,) = atm.provenance(("acabf",))
    assert row["version"] == "v2" and row["newer"] == ["v3"]
    # ...and the log says what it passed over
    (line,) = describe_forcing_provenance(atm, variables={"atmosphere": ("acabf",)})
    assert line.endswith("v2  NEWER ON DISK, NOT READ: v3")


def test_an_absent_tree_is_a_statement_not_a_silence(tmp_path):
    from icepack2_tools.forcing import ISMIP7Ocean, describe_forcing_provenance
    (line,) = describe_forcing_provenance(ISMIP7Ocean(data_root=str(tmp_path), scenario="ssp585"))
    assert line.endswith("ocean CESM2-WACCM ssp585: nothing on disk")


# --- guards on what the readers will and will not serve --------------------

def _write_ocean_chunk(vdir, esm, scenario, first, last, value, calendar=None):
    import cftime
    years = range(first, last + 1)
    time = ([cftime.DatetimeNoLeap(y, 1, 1) for y in years] if calendar == "noleap"
            else np.array([np.datetime64(f"{y}-01-01") for y in years]))
    ds = xr.Dataset(
        {"tf": (("time", "z", "y", "x"),
                np.stack([np.full((2, 3, 3), value + (y - first), dtype="float32") for y in years]))},
        coords={"time": time, "z": np.array([-30.0, -90.0]),
                "y": np.array([0.0, 8000.0, 16000.0]), "x": np.array([0.0, 8000.0, 16000.0])},
    )
    ds.to_netcdf(os.path.join(vdir, f"tf_AIS_{esm}_{scenario}_ocean_v3_{first}-{last}.nc"))


@pytest.fixture
def ocean_tree(tmp_path):
    esm, scenario = "CESM2-WACCM", "ssp585"
    vdir = tmp_path / esm / scenario / "ocean" / "tf" / "v3"
    vdir.mkdir(parents=True)
    _write_ocean_chunk(str(vdir), esm, scenario, 2280, 2289, 1.0)
    # 2290-2294 is the missing chunk; the last one is on a noleap calendar
    _write_ocean_chunk(str(vdir), esm, scenario, 2295, 2299, 5.0, calendar="noleap")
    return tmp_path


def test_the_ocean_holds_one_year_past_its_end_and_says_so(ocean_tree, capsys):
    from icepack2_tools.forcing import ISMIP7Ocean
    ocean = ISMIP7Ocean(data_root=str(ocean_tree))
    assert ocean.coverage("tf") == (2280, 2299)
    at = (np.array([8000.0]), np.array([8000.0]))
    assert ocean.get_thermal_forcing(2299, *at)[0] == pytest.approx(9.0)
    assert ocean.get_thermal_forcing(2300, *at)[0] == pytest.approx(9.0)
    ocean.get_thermal_forcing(2300, *at)
    assert capsys.readouterr().out.count("holding 2299, the last year on disk") == 1


@pytest.mark.parametrize("year, where", [(2279, "precedes the series"),
                                         (2292, "falls between the chunk files"),
                                         (2301, "past the end of the series")])
def test_the_ocean_refuses_a_year_it_does_not_hold(ocean_tree, year, where):
    r"""It used to serve the nearest chunk without a word, so a historical
    run starting before its ocean ran on the wrong decade and succeeded."""
    from icepack2_tools.forcing import ISMIP7Ocean
    with pytest.raises(FileNotFoundError, match=where):
        ISMIP7Ocean(data_root=str(ocean_tree)).get_thermal_forcing(year, np.array([0.0]), np.array([0.0]))


def test_a_file_with_no_time_slices_is_named(mri_tree):
    r"""The 2300 gradient files on the share were empty until May 2026 (#8)."""
    vdir = mri_tree / "MRI-ESM2-0" / "ssp585" / "GEMB-SDBN1-8000m" / "acabf" / "v1"
    empty = vdir / "acabf_AIS_MRI-ESM2-0_ssp585_GEMB-SDBN1-8000m_v1_2300.nc"
    xr.Dataset({"acabf": (("time", "y", "x"), np.zeros((0, 3, 3), dtype="float32"))},
               coords={"time": np.array([], dtype="datetime64[ns]"),
                       "y": np.arange(3.0), "x": np.arange(3.0)}).to_netcdf(empty)
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    with pytest.raises(ValueError, match=r"_2300\.nc has an empty time axis.*discussion #8"):
        atm._load_year("acabf", 2300)


def test_a_misfiled_year_is_not_an_available_year(mri_tree):
    r"""#41: ssp585 anomalies sat on the share named ``historical``. The
    loader builds the full name and would not find one, so neither may the
    availability gate that runs before it."""
    vdir = mri_tree / "MRI-ESM2-0" / "ssp585" / "GEMB-SDBN1-8000m" / "acabf" / "v1"
    (vdir / "acabf_AIS_MRI-ESM2-0_historical_GEMB-SDBN1-8000m_v1_2297.nc").write_bytes(b"")
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    assert atm.available_years("acabf") == [2296, 2298, 2299]


def test_the_zhou_climatology_is_found_at_whatever_version_each_variable_is(tmp_path):
    r"""Read on the cluster in September 2026: ``tf`` at v3, ``so`` and
    ``thetao`` at v4, the v4 being the fix for the July fault that shipped tf
    inside the so file. The path builder said v3 for all three."""
    from icepack2_tools.forcing import _oi_climatology_path
    base = tmp_path / "obs" / "ocean" / "climatology" / "zhou_annual_06_nov"
    (base / "tf" / "v3").mkdir(parents=True)
    (base / "so" / "v3").mkdir(parents=True)
    (base / "so" / "v4").mkdir(parents=True)
    tail = "_AIS_obs_ocean_climatology_zhou_annual_06_nov_"
    assert _oi_climatology_path(str(tmp_path), "tf", "06_nov").endswith(f"tf/v3/tf{tail}v3_1972-2024.nc")
    assert _oi_climatology_path(str(tmp_path), "so", "06_nov").endswith(f"so/v4/so{tail}v4_1972-2024.nc")
    # nothing on disk: the same missing path the callers already report
    assert _oi_climatology_path(str(tmp_path), "thetao", "06_nov").endswith(f"thetao/v3/thetao{tail}v3_1972-2024.nc")
    assert _oi_climatology_path(str(tmp_path), "tf", "30_sep").endswith("meltMIP/OI_Climatology_ismip8km_60m_tf_extrap.nc")


# --- the decode of a late year is silent ------------------------------------

def _serialization_warnings(call):
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = call()
    return result, [w for w in caught if issubclass(w.category, xr.SerializationWarning)]


def test_the_readers_decode_years_past_2262_without_a_serialization_warning(mri_tree, ocean_tree):
    r"""xarray decodes a standard-calendar axis to ``datetime64[ns]`` while
    the dates fit and falls back to ``cftime`` past 2262, warning on every
    open. The readers take only the year of a slice, so they ask for
    ``cftime`` up front: a 2015-2300 projection then reads its late chunks
    without printing a warning per chunk, and the values are the same."""
    atm = ISMIP7Atmosphere(data_root=str(mri_tree), esm="MRI-ESM2-0", scenario="ssp585", version="v1")
    field, warned = _serialization_warnings(lambda: atm._load_year("acabf", 2299))
    assert not warned, [str(w.message) for w in warned]
    assert np.allclose(np.asarray(field), 2.0)

    from icepack2_tools.forcing import ISMIP7Ocean
    ocean = ISMIP7Ocean(data_root=str(ocean_tree))
    at = (np.array([8000.0]), np.array([8000.0]))
    # 2289 sits in the standard-calendar chunk, 2299 in the noleap one
    for year, expect in ((2289, 10.0), (2299, 9.0)):
        tf, warned = _serialization_warnings(lambda: ocean.get_thermal_forcing(year, *at))
        assert not warned, [str(w.message) for w in warned]
        assert tf[0] == pytest.approx(expect)
