r"""The fluidity prior mean from a depth-averaged temperature raster
(icepack2_tools.rheology_prior): rate_factor of the raster's temperature at
every node, with the raster's fill values and stray out-of-range pixels
replaced by their nearest valid neighbour before interpolation.

Serial, a 20 km x 10 km rectangle under a 1 km GeoTIFF written here.
"""
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

firedrake = pytest.importorskip("firedrake")
from firedrake import (                                         # noqa: E402
    Function,
    FunctionSpace,
    RectangleMesh,
    SpatialCoordinate,
    VectorFunctionSpace,
)
from icepack.models.viscosity import rate_factor                # noqa: E402

from icepack2_tools.rheology_prior import (                     # noqa: E402
    T_MAX,
    T_MIN,
    clean_temperature,
    fluidity_prior_from_temperature_raster,
    temperature_on_space,
)

DX = 1000.0
X0, X1, Y0, Y1 = -2000.0, 22000.0, -2000.0, 12000.0
BAD = (10500.0, 5500.0)         # the pixel centre that gets the garbage


def _temperature(x, y):
    return 250.0 + 5e-4 * x + 2e-4 * y


def _raster(path, garbage):
    nx, ny = int((X1 - X0) / DX), int((Y1 - Y0) / DX)
    xc = X0 + DX * (np.arange(nx) + 0.5)
    yc = Y1 - DX * (np.arange(ny) + 0.5)          # north to south
    X, Y = np.meshgrid(xc, yc)
    T = _temperature(X, Y)
    if garbage is not None:
        i = int(np.argmin(np.abs(xc - BAD[0])))
        j = int(np.argmin(np.abs(yc - BAD[1])))
        T[j, i] = garbage
    with rasterio.open(
        path, "w", driver="GTiff", height=ny, width=nx, count=1,
        dtype="float64", crs="EPSG:3031", nodata=np.nan,
        transform=from_origin(X0, Y1, DX, DX),
    ) as dst:
        dst.write(T, 1)
    return str(path)


@pytest.fixture
def space():
    mesh = RectangleMesh(40, 20, 20e3, 10e3)
    return FunctionSpace(mesh, "CG", 1)


def _coords(Q):
    mesh = Q.mesh()
    return Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro


def test_clean_temperature_fills_from_the_nearest_valid_pixel():
    T = np.full((5, 5), 250.0)
    T[2, 2] = np.nan
    T[0, 0] = 4000.0
    T[4, 4] = -2400.0
    T[1, 3] = 260.0
    out, replaced = clean_temperature(T)
    assert replaced.sum() == 3
    assert np.all(np.isfinite(out))
    assert np.all((out >= T_MIN) & (out <= T_MAX))
    assert out[1, 3] == 260.0                     # valid pixels untouched
    assert out[2, 2] in (250.0, 260.0)            # a neighbour's value


def test_the_prior_is_the_rate_factor_of_the_raster_temperature(tmp_path, space):
    fn = _raster(tmp_path / "T.tif", garbage=None)
    A, info = fluidity_prior_from_temperature_raster(fn, space)
    xy = _coords(space)
    expected = rate_factor(_temperature(xy[:, 0], xy[:, 1]))
    # a linear field is reproduced exactly by bilinear interpolation
    assert np.allclose(A.dat.data_ro, expected, rtol=1e-9, atol=0.0)
    assert info["pixels_replaced"] == 0
    assert info["nodes_from_replaced"] == 0
    assert np.all(A.dat.data_ro > 0.0)


@pytest.mark.parametrize("garbage", [np.nan, 4481.7, -2409.7])
def test_a_bad_pixel_is_filled_and_the_nodes_it_reaches_are_counted(
        tmp_path, space, garbage):
    fn = _raster(tmp_path / "T.tif", garbage=garbage)
    T, info = temperature_on_space(fn, space)
    xy = _coords(space)
    exact = _temperature(xy[:, 0], xy[:, 1])
    assert info["pixels_replaced"] == 1
    assert 0 < info["nodes_from_replaced"] < 8
    assert np.all(np.isfinite(T.dat.data_ro))
    assert np.all((T.dat.data_ro >= T_MIN) & (T.dat.data_ro <= T_MAX))
    far = np.hypot(xy[:, 0] - BAD[0], xy[:, 1] - BAD[1]) > 1.5 * DX
    assert np.allclose(T.dat.data_ro[far], exact[far], rtol=1e-9, atol=0.0)
    # the filled pixel is a neighbour's temperature, so the nodes it reaches
    # are off by at most one pixel's worth of the gradient
    near = ~far
    assert np.all(np.abs(T.dat.data_ro[near] - exact[near]) <= 5e-4 * DX + 2e-4 * DX)
