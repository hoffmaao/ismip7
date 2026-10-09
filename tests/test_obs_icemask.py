"""One year of the Greene ice mask, cut from the MIPkit into a cached GeoTIFF
and classified against BedMachine on the same pixels (issue #167)."""

import datetime
import os

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")
rasterio = pytest.importorskip("rasterio")

from icepack2_tools import obs_icemask as OI  # noqa: E402


def _days(*dates):
    return [(datetime.date(*d) - datetime.date(1900, 1, 1)).days for d in dates]


def _grid(path, x, y, fields):
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("x", len(x))
        d.createDimension("y", len(y))
        for name, vals in (("x", x), ("y", y)):
            v = d.createVariable(name, "f8", (name,))
            v[:] = vals
            v.standard_name = f"projection_{name}_coordinate"
            v.units = "m"
        for name, (typ, vals) in fields.items():
            v = d.createVariable(name, typ, ("y", "x"))
            v[:] = vals
    return str(path)


def _kit(path):
    x = np.arange(0.0, 20000.0 + 1.0, 500.0)
    y = x[::-1].copy()
    X, _ = np.meshgrid(x, y)
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("x", len(x))
        d.createDimension("y", len(y))
        d.createDimension("greene_mask_time", 2)
        for name, vals in (("x", x), ("y", y)):
            v = d.createVariable(name, "i4", (name,))
            v[:] = vals
        t = d.createVariable("greene_mask_time", "i4", ("greene_mask_time",))
        t[:] = _days((2014, 3, 15), (2015, 3, 15))
        m = d.createVariable("icemask_greene", "i1", ("greene_mask_time", "y", "x"))
        m[0] = np.zeros(X.shape, "i1")
        m[1] = (X < 10000.0).astype("i1")
    return str(path)


def _bedmachine(path, shift=0.0):
    x = np.arange(-2000.0, 22000.0 + 1.0, 500.0) + shift
    y = x[::-1].copy()
    X, Y = np.meshgrid(x, y)
    mask = np.where(X < 8000.0, 2, np.where(X < 12000.0, 3, 0))
    mask = np.where((X >= 12000.0) & (Y > 15000.0), 1, mask)
    bed = np.where(X < 8000.0, 100.0, -500.0)
    bed = np.where(mask == 1, 50.0, bed)
    return _grid(path, x, y, {"mask": ("i1", mask), "bed": ("f4", bed),
                              "thickness": ("f4", np.where(X < 12000.0, 300.0, 0.0))})


@pytest.fixture
def kit_env(tmp_path, monkeypatch):
    kit = _kit(tmp_path / "AntarcticaObsISMIP7-v1.2.nc")
    monkeypatch.setenv("ISMIP7_OBS_KIT", kit)
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", str(tmp_path / "obs"))
    return tmp_path


def test_epoch_index_takes_the_march_mask_of_a_year_with_two():
    days = _days((1997, 10, 1), (2000, 3, 15), (2000, 10, 1), (2015, 3, 15))
    assert OI.epoch_index(days, 2000)[0] == 1
    assert OI.epoch_index(days, 2015) == (3, datetime.date(2015, 3, 15))
    with pytest.raises(ValueError):
        OI.epoch_index(days, 2016)


def test_the_year_is_cut_into_a_cached_tif(kit_env):
    tif = OI.icemask_tif(2015)
    assert tif == os.path.join(str(kit_env / "obs"), "icemask_cache",
                               "icemask_greene_2015_AntarcticaObsISMIP7-v1.2.tif")
    with rasterio.open(tif) as src:
        arr = src.read(1)
        assert src.tags()["date"] == "2015-03-15"
        left = src.bounds.left
    assert left == -250.0
    assert arr[:, :20].all() and not arr[:, 20:].any()
    mtime = os.path.getmtime(tif)
    assert OI.icemask_tif(2015) == tif and os.path.getmtime(tif) == mtime


def test_a_site_without_the_kit_reads_the_staged_tif(kit_env, monkeypatch):
    tif = OI.icemask_tif(2015)
    monkeypatch.delenv("ISMIP7_OBS_KIT")
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(kit_env / "nothing"))
    assert OI.icemask_source(2015) == (None, tif)
    with pytest.raises(FileNotFoundError):
        OI.icemask_source(2014)


def test_classes_on_the_mask_grid(kit_env):
    tif = OI.icemask_tif(2015)
    bm = _bedmachine(kit_env / "bm.nc")
    classes, transform = OI.read_classes(tif, bm)
    assert classes.shape == (41, 41)
    x = transform.c + 500.0 * (np.arange(41) + 0.5)
    y = transform.f - 500.0 * (np.arange(41) + 0.5)
    X, Y = np.meshgrid(x, y)
    expect = np.where(X < 10000.0, OI.ICE, OI.MARINE)
    expect = np.where((X >= 12000.0) & (Y > 15000.0), OI.LAND, expect)
    np.testing.assert_array_equal(classes, expect)


def test_bedmachine_off_the_mask_grid_is_refused(kit_env):
    tif = OI.icemask_tif(2015)
    bm = _bedmachine(kit_env / "bm_shifted.nc", shift=250.0)
    with pytest.raises(ValueError, match="pixel grid"):
        OI.read_classes(tif, bm)


def test_bedmachines_own_classes(tmp_path):
    bm = _bedmachine(tmp_path / "bm.nc")
    classes, transform = OI.bedmachine_classes(bm)
    assert classes.shape == (49, 49)
    x = transform.c + 500.0 * (np.arange(49) + 0.5)
    y = transform.f - 500.0 * (np.arange(49) + 0.5)
    X, Y = np.meshgrid(x, y)
    expect = np.where(X < 12000.0, OI.ICE, OI.MARINE)
    expect = np.where((X >= 12000.0) & (Y > 15000.0), OI.LAND, expect)
    np.testing.assert_array_equal(classes, expect)

