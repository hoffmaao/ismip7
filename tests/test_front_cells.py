r"""Front cells from BedMachine's own mask (issue #167).

Vertex sampling gives a cell on the ice front the mean of its vertex samples,
and a vertex on the front samples a blend of ice and water, so the front cells
hold a fraction of the front's thickness even when mesh edges follow the
front. Under ``vertex_front`` sampling BedMachine's mask decides which cells
hold ice at the marine front, and the front cells take its own thickness and
bed over their ice.

Serial, a 20 km square of 1 km cells on a BedMachine-shaped NetCDF: grounded
ice on a bed above sea level for x < 6 km, a 200 m floating shelf on a bed at
-500 m to the front, open ocean beyond.
"""

# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
import icepack2_tools.dual_friction  # noqa: F401,E402

import numpy as np                                              # noqa: E402
import netCDF4                                                  # noqa: E402
import pytest                                                   # noqa: E402
import rasterio                                                 # noqa: E402

from firedrake import (                                         # noqa: E402
    Function,
    FunctionSpace,
    RectangleMesh,
    SpatialCoordinate,
    VectorFunctionSpace,
)

from icepack2_tools import geometry as G                        # noqa: E402

SHELF_H = 200.0
GROUNDED_H = 800.0
PX = 500.0
X = np.arange(-2000.0, 22000.0 + 1.0, PX)                        # pixel centres
Y = X[::-1].copy()


def _bedmachine(path, front_x=12000.0, nunatak=False):
    XX, YY = np.meshgrid(X, Y)
    grounded = XX < 6000.0
    ice = XX < front_x
    mask = np.where(grounded, 2, np.where(ice, 3, 0))
    thk = np.where(grounded, GROUNDED_H, np.where(ice, SHELF_H, 0.0))
    bed = np.where(grounded, 100.0, -500.0)
    if nunatak:
        rock = (np.abs(XX - 3000.0) < 1000.0) & (np.abs(YY - 10000.0) < 1000.0)
        mask = np.where(rock, 1, mask)
        thk = np.where(rock, 0.0, thk)
        bed = np.where(rock, 900.0, bed)
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("x", len(X))
        d.createDimension("y", len(Y))
        for name, vals in (("x", X), ("y", Y)):
            v = d.createVariable(name, "f8", (name,))
            v[:] = vals
            v.standard_name = f"projection_{name}_coordinate"
            v.units = "m"
        for name, typ, vals in (("mask", "i1", mask), ("thickness", "f4", thk),
                                ("bed", "f4", bed),
                                ("surface", "f4", np.where(grounded, bed + thk, 0.1 * thk))):
            v = d.createVariable(name, typ, ("y", "x"))
            v[:] = vals.astype(typ)
    return str(path)


def _spaces(offset=0.0):
    mesh = RectangleMesh(20, 20, 20e3, 20e3)
    if offset:
        mesh.coordinates.dat.data[:, 0] += offset
    Q = FunctionSpace(mesh, "CG", 1)
    Q0 = FunctionSpace(mesh, "DG", 0)
    xc = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro
    return Q, Q0, xc


def _vertex(bm, Q0, Q):
    b, H, counts = G.sample_bed_thickness(bm, Q0, Q, method="vertex")
    assert counts is None
    return b.dat.data_ro.copy(), H.dat.data_ro.copy()


def test_vertex_sampling_thins_the_front_cells(tmp_path):
    bm = _bedmachine(tmp_path / "bm.nc")
    Q, Q0, xc = _spaces()
    _, H = _vertex(bm, Q0, Q)
    front = (xc[:, 0] > 11e3) & (xc[:, 0] < 12e3)
    # a vertex on the front does not sample the full shelf, so every front
    # cell is thinner than the shelf
    assert H[front].max() < SHELF_H - 1.0


def test_front_cells_take_the_shelf_thickness(tmp_path):
    bm = _bedmachine(tmp_path / "bm.nc")
    Q, Q0, xc = _spaces()
    b0, H0 = _vertex(bm, Q0, Q)
    b, H, _ = G.sample_bed_thickness(bm, Q0, Q, method="vertex")
    counts = G.front_cells(H, b, bm)
    H, b = H.dat.data_ro, b.dat.data_ro
    front = (xc[:, 0] > 11e3) & (xc[:, 0] < 12e3)
    water = xc[:, 0] > 12e3
    interior = xc[:, 0] < 10e3
    np.testing.assert_allclose(H[front], SHELF_H)
    np.testing.assert_allclose(b[front], -500.0)
    assert np.all(H[water] == 0.0)
    np.testing.assert_array_equal(H[interior], H0[interior])
    np.testing.assert_array_equal(b[interior], b0[interior])
    assert counts["rebuilt"] == front.sum() and counts["water"] == water.sum()


def test_a_cell_holds_ice_when_most_of_it_is_ice(tmp_path):
    # Pixels are centred on multiples of 500 m, so the ice ends at the pixel
    # edge at 12.25 km. With the mesh moved 50 m east the 12.05-13.05 km
    # column holds 20 % ice, so both its triangles are water, and the column
    # before it is the front.
    bm = _bedmachine(tmp_path / "bm.nc", front_x=12300.0)
    Q, Q0, xc = _spaces(offset=50.0)
    b, H, _ = G.sample_bed_thickness(bm, Q0, Q, method="vertex")
    G.front_cells(H, b, bm)
    H = H.dat.data_ro
    assert np.all(H[xc[:, 0] > 12.05e3] == 0.0)
    np.testing.assert_allclose(H[(xc[:, 0] > 11.05e3) & (xc[:, 0] < 12.05e3)], SHELF_H)
    # moved 550 m west, the 11.45-12.45 km column holds 80 % ice: it holds the
    # shelf's full thickness, and holds water itself, so it is rebuilt
    Q, Q0, xc = _spaces(offset=-550.0)
    b, H, _ = G.sample_bed_thickness(bm, Q0, Q, method="vertex")
    G.front_cells(H, b, bm)
    H = H.dat.data_ro
    np.testing.assert_allclose(H[(xc[:, 0] > 11.45e3) & (xc[:, 0] < 12.45e3)], SHELF_H)
    assert np.all(H[xc[:, 0] > 12.45e3] == 0.0)


def test_a_nunatak_is_left_alone(tmp_path):
    bm = _bedmachine(tmp_path / "bm.nc", nunatak=True)
    Q, Q0, xc = _spaces()
    b0, H0 = _vertex(bm, Q0, Q)
    b, H, _ = G.sample_bed_thickness(bm, Q0, Q, method="vertex")
    G.front_cells(H, b, bm)
    near = (np.abs(xc[:, 0] - 3e3) < 2e3) & (np.abs(xc[:, 1] - 10e3) < 2e3)
    np.testing.assert_array_equal(H.dat.data_ro[near], H0[near])
    np.testing.assert_array_equal(b.dat.data_ro[near], b0[near])


def test_vertex_front_sampling_rebuilds_the_front(tmp_path):
    bm = _bedmachine(tmp_path / "bm.nc")
    Q, Q0, xc = _spaces()
    b, H, counts = G.sample_bed_thickness(bm, Q0, Q, method="vertex_front")
    front = (xc[:, 0] > 11e3) & (xc[:, 0] < 12e3)
    assert counts["rebuilt"] == front.sum()
    np.testing.assert_allclose(H.dat.data_ro[front], SHELF_H)
    with pytest.raises(ValueError):
        G.sample_bed_thickness(bm, Q, Q, method="vertex_front")


def test_the_lattice_mean_of_a_ramp_is_the_cell_centroid(tmp_path):
    path = tmp_path / "ramp.nc"
    XX, _ = np.meshgrid(X, Y)
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("x", len(X))
        d.createDimension("y", len(Y))
        for name, vals in (("x", X), ("y", Y)):
            v = d.createVariable(name, "f8", (name,))
            v[:] = vals
            v.standard_name = f"projection_{name}_coordinate"
        v = d.createVariable("ramp", "f4", ("y", "x"))
        v[:] = XX.astype("f4")
    Q, Q0, xc = _spaces()
    mean = G.raster_cell_mean(rasterio.open(f"netcdf:{path}:ramp"), Q0)
    assert np.abs(mean.dat.data_ro - xc[:, 0]).max() < PX / 2
