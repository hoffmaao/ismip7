r"""The ISMIP7 ice-shelf collapse mask reaches the transport's cells.

Protocol path C (discussions #30, #33): the mask applies to floating ice
only, exists for the SSP scenarios and not for historical or OCX. This
checks the reader finds a versioned mask file and samples it by nearest
pixel, and that the forcing callback fills ``ctx["collapse"]`` on the
geometry cells when the run provides it. The removal itself is the
transport's (``ISMIP7_FRACTURE=mask``), exercised by the forward.
"""
import os

import numpy as np
import pytest

xr = pytest.importorskip("xarray")

from icepack2_tools.forcing import ISMIP7Fracture, make_forcing_callback  # noqa: E402


@pytest.fixture
def mask_tree(tmp_path):
    d = tmp_path / "CESM2-WACCM" / "ssp585" / "fracture" / "v2.1"
    d.mkdir(parents=True)
    x = np.array([0.0, 8000.0, 16000.0]); y = np.array([0.0, 8000.0])
    m = np.zeros((2, 2, 3), dtype="int8")        # (time, y, x)
    m[0, 1, 2] = 1                                # 2050: the cell at (16 km, 8 km)
    m[1, :, :] = 1                                # 2100: everything
    ds = xr.Dataset({"ice_shelf_collapse_mask": (("time", "y", "x"), m)},
                    coords={"time": np.array([2050, 2100]), "y": y, "x": x})
    ds.to_netcdf(d / "ice_shelf_collapse_mask_cesm2waccm_ssp585_ismip7_8km-v2.1.nc")
    return tmp_path


def test_reader_samples_the_versioned_mask(mask_tree):
    fr = ISMIP7Fracture(data_root=str(mask_tree), esm="CESM2-WACCM", scenario="ssp585").load()
    # points inside the pixel-centre extent: nearest pixel wins; a point
    # beyond the outermost centres is outside the grid and reads as 0
    xs = np.array([100.0, 15900.0, 8100.0, 20000.0]); ys = np.array([100.0, 7900.0, 100.0, 100.0])
    assert fr.get_collapse_mask(2050, xs, ys).tolist() == [0.0, 1.0, 0.0, 0.0]
    assert fr.get_collapse_mask(2100, xs, ys).tolist() == [1.0, 1.0, 1.0, 0.0]
    assert fr.get_collapse_mask(2030, xs, ys).tolist() == [0.0, 1.0, 0.0, 0.0]   # nearest year


def test_the_year_is_read_however_the_file_encodes_time(tmp_path):
    r"""A plain-year axis is only one way to write the file. With CF units
    xarray decodes the axis to dates, and ``sel(time=2050)`` raises on it."""
    d = tmp_path / "CESM2-WACCM" / "ssp585" / "fracture" / "v2.1"
    d.mkdir(parents=True)
    m = np.zeros((2, 1, 1), dtype="int8")
    m[1] = 1
    ds = xr.Dataset({"ice_shelf_collapse_mask": (("time", "y", "x"), m)},
                    coords={"time": np.array(["2050-07-01", "2100-07-01"], dtype="datetime64[ns]"),
                            "y": np.array([0.0]), "x": np.array([0.0])})
    ds.to_netcdf(d / "ice_shelf_collapse_mask_cesm2waccm_ssp585_ismip7_8km-v2.1.nc")
    fr = ISMIP7Fracture(data_root=str(tmp_path), esm="CESM2-WACCM", scenario="ssp585").load()
    at = (np.array([0.0]), np.array([0.0]))
    assert fr.get_collapse_mask(2060, *at).tolist() == [0.0]
    assert fr.get_collapse_mask(2090, *at).tolist() == [1.0]
    assert fr.provenance()[0]["version"] == "v2.1"


def test_the_mask_is_found_by_name_among_the_files_other_variables(tmp_path):
    r"""The real files (read on the cluster, September 2026) hold ``mask``
    beside a scalar ``mapping`` and 2-D ``lon``/``lat``. Taking the first data
    variable worked only because ``mask`` happens to come first."""
    d = tmp_path / "CESM2-WACCM" / "ssp585" / "fracture" / "v2.1"
    d.mkdir(parents=True)
    lon = np.full((1, 2), 170.0)
    ds = xr.Dataset({"mapping": ((), np.int32(0)),
                     "lon": (("y", "x"), lon),
                     "mask": (("time", "y", "x"), np.array([[[0, 1]]], dtype="int8"),
                              {"standard_name": "ice_shelf_collapse_mask", "units": "1"})},
                    coords={"time": ("time", np.array([2050], dtype="int32"), {"units": "year"}),
                            "y": np.array([0.0]), "x": np.array([0.0, 8000.0])})
    ds.to_netcdf(d / "ice_shelf_collapse_mask_cesm2waccm_ssp585_ismip7_8km-v2.1.nc")
    fr = ISMIP7Fracture(data_root=str(tmp_path), esm="CESM2-WACCM", scenario="ssp585").load()
    assert fr.get_collapse_mask(2299, np.array([0.0, 8000.0]), np.array([0.0, 0.0])).tolist() == [0.0, 1.0]


def test_callback_fills_the_collapse_cells(mask_tree):
    fr = ISMIP7Fracture(data_root=str(mask_tree), esm="CESM2-WACCM", scenario="ssp585").load()
    cb = make_forcing_callback(fracture=fr)
    ctx = {"geom_xy": (np.array([100.0, 15900.0]), np.array([100.0, 7900.0])),
           "collapse": np.zeros(2, dtype=bool)}
    cb(ctx, 2050.0)
    assert list(ctx["collapse"]) == [False, True]
    ctx.pop("collapse")
    cb(ctx, 2050.0)                                               # a run without the knob: nothing to fill, no error


def test_missing_mask_means_no_collapse(tmp_path):
    fr = ISMIP7Fracture(data_root=str(tmp_path), esm="CESM2-WACCM", scenario="historical").load()
    assert not fr.get_collapse_mask(2000, np.zeros(3), np.zeros(3)).any()


def _mask(root, esm, scenario, version):
    r"""An empty collapse mask named as the focus group names them."""
    tag = {"CESM2-WACCM": "cesm2waccm", "MRI-ESM2-0": "mriesm20"}[esm]
    d = root / esm / scenario / "fracture" / version
    d.mkdir(parents=True, exist_ok=True)
    ds = xr.Dataset({"mask": (("time", "y", "x"), np.zeros((1, 2, 2), dtype="int8"))},
                    coords={"time": np.array([2100]), "y": np.array([0.0, 8000.0]),
                            "x": np.array([0.0, 8000.0])})
    ds.to_netcdf(d / f"ice_shelf_collapse_mask_{tag}_{scenario}_ismip7_8km-{version}.nc")


def test_a_mask_mode_run_refuses_a_replaced_mask(tmp_path):
    r"""The focus group replaced the MRI-ESM2-0 ssp585 v1 mask with v2 on Globus
    on 22 September 2026 while the mirror still served v1 (discussion #30). A
    tree synced from the mirror alone has to stop a mask-mode run, and the
    message has to say where the fix is; a tree holding v2 runs."""
    _mask(tmp_path, "MRI-ESM2-0", "ssp585", "v1")
    fr = ISMIP7Fracture(data_root=str(tmp_path), esm="MRI-ESM2-0", scenario="ssp585").load()
    assert fr.version() == "v1"
    with pytest.raises(RuntimeError, match=r"needs v2 or newer.*on Globus.*"
                                           r"download_forcing\.py --scenarios --esm MRI-ESM2-0 "
                                           r"--scenario ssp585"):
        fr.check_min_version()
    _mask(tmp_path, "MRI-ESM2-0", "ssp585", "v2")
    fr = ISMIP7Fracture(data_root=str(tmp_path), esm="MRI-ESM2-0", scenario="ssp585").load()
    assert fr.version() == "v2"
    fr.check_min_version()


def test_the_minimum_is_per_esm_and_scenario(mask_tree):
    r"""MRI-ESM2-0 ssp126 has no fix, so its v1 is current; CESM2-WACCM is at
    v2.1 everywhere, and its v2 is the replaced one."""
    _mask(mask_tree, "MRI-ESM2-0", "ssp126", "v1")
    ISMIP7Fracture(data_root=str(mask_tree), esm="MRI-ESM2-0", scenario="ssp126").load().check_min_version()
    ISMIP7Fracture(data_root=str(mask_tree), esm="CESM2-WACCM", scenario="ssp585").load().check_min_version()
    _mask(mask_tree, "CESM2-WACCM", "ssp370", "v2")
    with pytest.raises(RuntimeError, match=r"needs v2\.1 or newer"):
        ISMIP7Fracture(data_root=str(mask_tree), esm="CESM2-WACCM", scenario="ssp370").load().check_min_version()
