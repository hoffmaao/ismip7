r"""Core 11 reads the ISMIP7 OCX product, and says so when it does not.

The observation-constrained experiment has no ESM and no scenario, and the
focus groups decided its tree need only be consistent with itself (discussion
#41, item 6): directories run ``OCX/<source>/<product>/...``, scenario first,
while filenames keep the source before ``OCX``. Its Antarctic ocean cites no
source at all: four expert-judgment scenarios, each one directory holding
``so``, ``tf`` and ``thetao`` side by side, versioned ``v1`` since 28 August
2026. The readers could address none of it, so core 11 ran on a RACMO2.4p1
and OI-climatology stopgap with the real product on disk.

Discussion #48 is why reading it needs a tripwire: the OCX ``main`` thermal
forcing disagrees with the climatology around Mertz, and every K here is
fitted to the climatology.
"""
import os
import sys

import numpy as np
import pytest

xr = pytest.importorskip("xarray")

from icepack2_tools.forcing import (OCX, OCX_ATMOSPHERE_SOURCE, ISMIP7Atmosphere,  # noqa: E402
                                    ISMIP7Ocean, describe_forcing_provenance)
from icepack2_tools.runconfig import ocx_forcing, ocx_ocean  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

XY = {"y": np.array([0.0, 8000.0, 16000.0]), "x": np.array([0.0, 8000.0, 16000.0])}


@pytest.fixture
def ocx_tree(tmp_path):
    adir = tmp_path / "OCX" / OCX_ATMOSPHERE_SOURCE / "SDBN1-8000m" / "acabf" / "v1"
    adir.mkdir(parents=True)
    for year, value in ((1979, 1.0), (1980, 2.0)):
        times = np.array([np.datetime64(f"{year}-{m:02d}-15") for m in range(1, 13)])
        xr.Dataset({"acabf": (("time", "y", "x"), np.full((12, 3, 3), value, dtype="float32"))},
                   coords={"time": times, **XY}).to_netcdf(
            adir / f"acabf_AIS_{OCX_ATMOSPHERE_SOURCE}_OCX_SDBN1-8000m_v1_{year}.nc")
    odir = tmp_path / "OCX" / "ocean" / "main" / "v1"
    odir.mkdir(parents=True)
    years = np.arange(1950, 2026)
    for var, base in (("tf", 1.0), ("so", 34.0), ("thetao", -1.0)):
        field = (base + (years - 1950)[:, None, None, None] * 0.01) * np.ones((1, 2, 3, 3))
        xr.Dataset({var: (("time", "z", "y", "x"), field.astype("float32"))},
                   coords={"time": np.array([np.datetime64(f"{y}-01-01") for y in years]),
                           "z": np.array([-30.0, -90.0]), **XY}).to_netcdf(
            odir / f"{var}_AIS_OCX_ocean_main_v1_1950-2025.nc")
    return tmp_path


def test_the_atmosphere_reader_finds_the_scenario_first_tree(ocx_tree):
    atm = ISMIP7Atmosphere(data_root=str(ocx_tree), esm=OCX_ATMOSPHERE_SOURCE, scenario=OCX)
    # pinned v2 is not there, so the highest on disk, v1, is what it opens
    assert atm.available_years("acabf") == [1979, 1980]
    assert np.allclose(np.asarray(atm._load_year("acabf", 1980)), 2.0)
    # the old driver asked for <ESM>/ocx/..., which no tree ever had
    assert ISMIP7Atmosphere(data_root=str(ocx_tree), scenario="ocx").available_years("acabf") == []


def test_the_ocean_reader_takes_one_variable_from_the_shared_directory(ocx_tree):
    ocean = ISMIP7Ocean(data_root=str(ocx_tree), scenario=OCX)
    assert ocean.coverage("tf") == ocean.coverage("so") == (1950, 2025)
    assert [os.path.basename(p) for _, _, p in ocean.spans("tf")] == ["tf_AIS_OCX_ocean_main_v1_1950-2025.nc"]
    at = (np.array([8000.0]), np.array([8000.0]))
    assert ocean.get_thermal_forcing(2000, *at)[0] == pytest.approx(1.5)
    assert ocean.get_salinity(1950, *at)[0] == pytest.approx(34.0)
    with pytest.raises(FileNotFoundError, match="precedes the series"):
        ocean.get_thermal_forcing(1949, *at)


def test_an_absent_variant_is_absent_and_a_misspelt_one_is_refused(ocx_tree):
    assert ISMIP7Ocean(data_root=str(ocx_tree), scenario=OCX, variant="warm").coverage("tf") is None
    with pytest.raises(ValueError, match="main"):
        ISMIP7Ocean(data_root=str(ocx_tree), scenario=OCX, variant="hot")


def test_the_report_says_which_ocx_product_a_run_opened(ocx_tree):
    atm = ISMIP7Atmosphere(data_root=str(ocx_tree), esm=OCX_ATMOSPHERE_SOURCE, scenario=OCX)
    ocean = ISMIP7Ocean(data_root=str(ocx_tree), scenario=OCX)
    lines = describe_forcing_provenance(atm, ocean, variables={"atmosphere": ("acabf",)})
    assert [ln.split(": ", 1)[1] for ln in lines] == [
        "atmosphere acabf RACMO2.3p2-ERA OCX SDBN1-8000m v1",
        "ocean tf expert-judgment OCX ocean/main v1",
        "ocean so expert-judgment OCX ocean/main v1",
    ]


def test_the_stopgap_has_to_be_asked_for(monkeypatch):
    monkeypatch.delenv("ISMIP7_OCX_FORCING", raising=False)
    monkeypatch.delenv("ISMIP7_OCX_OCEAN", raising=False)
    assert (ocx_forcing(), ocx_ocean()) == ("protocol", "main")
    monkeypatch.setenv("ISMIP7_OCX_FORCING", "stopgap")
    monkeypatch.setenv("ISMIP7_OCX_OCEAN", "cold")
    assert (ocx_forcing(), ocx_ocean()) == ("stopgap", "cold")
    monkeypatch.setenv("ISMIP7_OCX_FORCING", "racmo")
    with pytest.raises(ValueError, match="protocol"):
        ocx_forcing()


def test_the_driver_refuses_a_product_that_does_not_cover_the_run(ocx_tree, monkeypatch):
    pytest.importorskip("firedrake")
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(ocx_tree))
    monkeypatch.delenv("ISMIP7_OCX_OCEAN", raising=False)
    import importlib
    ocx = importlib.import_module("projections.ocx")
    atm, ocean = ocx.protocol_forcing(1979.0, 1981.0)          # years 1979, 1980
    assert (atm.esm, ocean.variant) == (OCX_ATMOSPHERE_SOURCE, "main")
    ocx.protocol_forcing(1979.0, 1982.0)                       # 1981 is the one held year
    with pytest.raises(FileNotFoundError, match="ISMIP7_OCX_FORCING=stopgap"):
        ocx.protocol_forcing(1979.0, 2026.0)
    monkeypatch.setenv("ISMIP7_OCX_OCEAN", "vary")
    with pytest.raises(FileNotFoundError, match=r"data/OCX/ocean/vary/"):
        ocx.protocol_forcing(1979.0, 1981.0)


# --- the discussion #48 tripwire -------------------------------------------

def test_a_region_whose_ocx_melt_is_off_is_flagged():
    pytest.importorskip("firedrake")
    pytest.importorskip("rasterio")
    import check_melt_bound as cmb
    x = np.array([10e3, 20e3, 300e3, 310e3, 900e3])
    y = np.array([10e3, 20e3, 10e3, 20e3, -300e3])
    ids, centres = cmb.block_ids(x, y)
    assert ids[0] == ids[1] != ids[2] == ids[3] != ids[4]
    assert centres[int(ids[0])] == (128e3, 128e3) and centres[int(ids[4])] == (896e3, -384e3)

    area = np.full(5, 64e6)                                    # 8 km pixels
    floating = np.array([True, True, True, True, False])
    clim = np.array([10.0, 10.0, 10.0, 10.0, 99.0])            # m/yr of ice
    ocx = np.array([10.0, 11.0, 5.0, 5.0, 0.0])                # the second block halves, as Mertz did
    ref = cmb.melt_by_group(ids, clim, area, floating)
    new = cmb.melt_by_group(ids, ocx, area, floating)
    assert set(ref) == {int(ids[0]), int(ids[2])}              # the grounded dof is in no total
    assert ref[int(ids[0])] == pytest.approx(2 * 10.0 * 64e6 * 917.0 / 1e12)
    assert cmb.off_by_more_than(ref, new, tol=0.25, floor_gt=0.1) == [int(ids[2])]
    # under the floor a ratio of two small totals is not worth reading
    assert cmb.off_by_more_than(ref, new, tol=0.25, floor_gt=5.0) == []


def test_the_check_names_the_ocx_version_it_was_measured_on(ocx_tree, monkeypatch, capsys):
    pytest.importorskip("firedrake")
    pytest.importorskip("rasterio")
    import argparse
    import check_melt_bound as cmb
    from icepack2_tools.forcing import OCEAN_VERSION
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(ocx_tree))
    # one basin everywhere: the IMBIE2 grid is not what this is about
    monkeypatch.setattr(cmb.cm, "_grid_interp",
                        lambda path, var, x, y, draft=None: np.ones(len(x)))
    n = 3
    g = {"x": XY["x"], "y": XY["y"], "draft": np.full(n, -60.0),
         "sin_a": np.full(n, 1e-2), "floating": np.ones(n, dtype=bool),
         "area": np.full(n, 64e6), "K": np.full(n, 1e-4),
         "tf": np.full(n, 1.0), "sal": np.full(n, 34.0)}          # the tree's 1950
    a = argparse.Namespace(ocx="main", ocx_years="1950", ocx_tol=0.25, ocx_floor=0.0)

    assert cmb.compare_with_ocx(a, g) == 0
    said = capsys.readouterr().out
    # The tree holds v1 alone, so the reader fell back from its pin, which is
    # the case at every site until upstream ships the pinned version.
    assert OCEAN_VERSION != "v1"
    assert "ocean tf expert-judgment OCX ocean/main v1" in said
    assert "ocean so expert-judgment OCX ocean/main v1" in said
