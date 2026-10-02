r"""Under BedMachine's subglacial-lake mask the bed is the lake floor, so
``b + H`` would sink the model surface by the water column. The correction
raises the bed to the ice base ``s - H`` over the lake and nowhere else.

Serial, a 20 km square mesh on a small BedMachine-shaped NetCDF with a 4 km
lake of 400 m of water under 4000 m of ice.
"""

# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
import icepack2_tools.dual_friction  # noqa: F401,E402

import numpy as np                                              # noqa: E402
import netCDF4                                                  # noqa: E402
import rasterio                                                 # noqa: E402
import pytest                                                   # noqa: E402

from firedrake import (                                         # noqa: E402
    Function,
    FunctionSpace,
    RectangleMesh,
    SpatialCoordinate,
    VectorFunctionSpace,
)

from icepack2_tools.geometry import (                           # noqa: E402
    raise_bed_to_lake_ice_base,
    sample_to_geometry,
)

LAKE = dict(bed=-900.0, thickness=4000.0, surface=3500.0, mask=4)
LAND = dict(bed=100.0, thickness=3000.0, surface=3100.0, mask=2)


def _bedmachine(path, with_lake=True):
    x = np.arange(-2000.0, 22000.0 + 1.0, 500.0)
    y = x[::-1].copy()
    X, Y = np.meshgrid(x, y)
    lake = with_lake & (X > 8000.0) & (X < 12000.0) & (Y > 8000.0) & (Y < 12000.0)
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("x", len(x))
        d.createDimension("y", len(y))
        for name, vals in (("x", x), ("y", y)):
            v = d.createVariable(name, "f8", (name,))
            v[:] = vals
            v.standard_name = f"projection_{name}_coordinate"
            v.units = "m"
        for name in ("bed", "thickness", "surface", "mask"):
            typ = "i1" if name == "mask" else "f4"
            v = d.createVariable(name, typ, ("y", "x"))
            v[:] = np.where(lake, LAKE[name], LAND[name]).astype(typ)
    return str(path)


def _geometry(bm_fn):
    mesh = RectangleMesh(20, 20, 20e3, 20e3)
    Q = FunctionSpace(mesh, "CG", 1)
    Q0 = FunctionSpace(mesh, "DG", 0)
    b = sample_to_geometry(rasterio.open(f"netcdf:{bm_fn}:bed"), Q0, Q)
    H = sample_to_geometry(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q0, Q)
    xc = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro
    return Q, Q0, b, H, xc


def test_the_bed_over_the_lake_becomes_the_ice_base(tmp_path):
    bm_fn = _bedmachine(tmp_path / "bm.nc")
    Q, Q0, b, H, xc = _geometry(bm_fn)
    b0 = b.dat.data_ro.copy()
    n = raise_bed_to_lake_ice_base(b, H, bm_fn, Q0, Q)
    inside = (np.abs(xc[:, 0] - 10e3) < 1e3) & (np.abs(xc[:, 1] - 10e3) < 1e3)
    far = (np.abs(xc[:, 0] - 10e3) > 4e3) | (np.abs(xc[:, 1] - 10e3) > 4e3)
    ice_base = LAKE["surface"] - LAKE["thickness"]
    assert np.allclose(b.dat.data_ro[inside], ice_base, atol=1e-3)
    assert np.allclose(b.dat.data_ro[inside] + H.dat.data_ro[inside], LAKE["surface"], atol=1e-3)
    assert np.array_equal(b.dat.data_ro[far], b0[far])
    changed = np.abs(b.dat.data_ro - b0) > 1e-6
    assert n == changed.sum() and inside.sum() <= n < 50


def test_no_lake_changes_nothing(tmp_path):
    bm_fn = _bedmachine(tmp_path / "bm.nc", with_lake=False)
    Q, Q0, b, H, xc = _geometry(bm_fn)
    b0 = b.dat.data_ro.copy()
    assert raise_bed_to_lake_ice_base(b, H, bm_fn, Q0, Q) == 0
    assert np.array_equal(b.dat.data_ro, b0)
