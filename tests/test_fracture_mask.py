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
