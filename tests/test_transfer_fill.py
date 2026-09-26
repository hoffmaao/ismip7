r"""Cross-mesh interpolation with a stated fill for the dofs the source mesh
does not cover (icepack2_tools.transfer.interpolate_with_fill).

The production case: a MAP inverted on a buffer-0 mesh transferred onto the
buffered production mesh. Every target dof in the 20 km ocean ring is
outside the source, and the fluidity prior must not become zero there.
"""
import numpy as np
import pytest

fd = pytest.importorskip("firedrake")

from icepack2_tools.transfer import (  # noqa: E402
    STRICT_TOLERANCE, interpolate_with_fill, strict_location,
)


def _meshes():
    # Target vertices at multiples of 0.3, so none sits exactly on the source
    # boundary x = 1 or y = 1 and the count is unambiguous.
    source = fd.UnitSquareMesh(4, 4)
    target = fd.RectangleMesh(5, 5, 1.5, 1.5)
    return source, target


def _outside(target):
    xy = target.coordinates.dat.data_ro
    return (xy[:, 0] > 1.0 + 1e-9) | (xy[:, 1] > 1.0 + 1e-9)


def test_a_scalar_fill_lands_on_every_dof_the_source_misses():
    source, target = _meshes()
    x, y = fd.SpatialCoordinate(source)
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(1 + x + 2 * y)
    g = fd.Function(fd.FunctionSpace(target, "CG", 1), name="g")
    n_missing, n_total, n_clamped = interpolate_with_fill(g, f, 7.0)
    outside = _outside(target)
    assert n_clamped == 0
    assert n_missing == int(outside.sum()) == 20
    assert n_total == 36
    vals = g.dat.data_ro
    assert np.all(vals[outside] == 7.0)
    xy = target.coordinates.dat.data_ro
    inside = ~outside
    assert np.allclose(vals[inside], 1 + xy[inside, 0] + 2 * xy[inside, 1], atol=1e-12)
    assert not np.isnan(vals).any()


def test_a_function_fill_supplies_the_uncovered_dofs_of_a_vector_field():
    source, target = _meshes()
    x, y = fd.SpatialCoordinate(source)
    V_s = fd.VectorFunctionSpace(source, "CG", 1)
    u = fd.Function(V_s).interpolate(fd.as_vector((x, y)))
    V_t = fd.VectorFunctionSpace(target, "CG", 1)
    X, Y = fd.SpatialCoordinate(target)
    raster = fd.Function(V_t).interpolate(fd.as_vector((X + 100.0, Y + 200.0)))
    w = fd.Function(V_t, name="velocity_obs")
    n_missing, n_total, n_clamped = interpolate_with_fill(w, u, raster)
    assert n_clamped == 0
    outside = _outside(target)
    assert n_missing == 20 and n_total == 36
    vals = w.dat.data_ro
    assert np.allclose(vals[outside], raster.dat.data_ro[outside])
    xy = target.coordinates.dat.data_ro
    assert np.allclose(vals[~outside], xy[~outside], atol=1e-12)


def test_a_fill_on_another_space_is_refused():
    source, target = _meshes()
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(fd.Constant(1.0))
    g = fd.Function(fd.FunctionSpace(target, "CG", 1))
    wrong = fd.Function(fd.FunctionSpace(target, "DG", 0))
    with pytest.raises(ValueError, match="target's function space"):
        interpolate_with_fill(g, f, wrong)


def test_a_same_mesh_transfer_misses_nothing():
    source, _ = _meshes()
    x, _y = fd.SpatialCoordinate(source)
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(x)
    g = fd.Function(fd.FunctionSpace(source, "CG", 1))
    n_missing, n_total, n_clamped = interpolate_with_fill(g, f, 7.0)
    assert n_missing == 0 and n_total == 25 and n_clamped == 0
    assert np.allclose(g.dat.data_ro, f.dat.data_ro)


def _shifted_meshes(shift=0.03):
    # The target reaches past the source on every side, and its first row of
    # dofs sits ``shift`` outside the outline: closer than half a source cell
    # (0.0625 in reference L1 distance), which Firedrake's default tolerance
    # accepts as inside, and further than the strict tolerance rejects.
    source = fd.UnitSquareMesh(8, 8)
    target = fd.RectangleMesh(12, 12, 1.5, 1.5)
    target.coordinates.dat.data[:] -= shift
    return source, target


def test_a_point_within_the_location_tolerance_is_a_fill_not_an_extrapolation():
    r"""Firedrake locates a target point up to half a reference cell outside
    a boundary cell (mesh.tolerance 0.5) and extrapolates that cell's linear
    basis there: the 2 km MAPs onto the 1 km mesh gave a fluidity prior of
    -218 from a source whose minimum was 1. The transfer locates strictly,
    so the band beyond the outline is a fill like the rest of the ring and
    nothing is left for the clamp."""
    source = fd.UnitSquareMesh(4, 4)
    x, _y = fd.SpatialCoordinate(source)
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(1 + 10 * x)
    target = fd.RectangleMesh(11, 11, 1.1, 1.1)
    g = fd.Function(fd.FunctionSpace(target, "CG", 1))
    n_missing, n_total, n_clamped = interpolate_with_fill(g, f, 7.0)
    xy = target.coordinates.dat.data_ro
    vals = g.dat.data_ro
    beyond = np.isclose(xy[:, 0], 1.1) & (xy[:, 1] <= 1.0 + 1e-9)
    assert beyond.sum() == 11
    assert np.all(vals[beyond] == 7.0)   # 12 by extrapolation, 11 by the clamp
    outside = (xy[:, 0] > 1.0 + 1e-9) | (xy[:, 1] > 1.0 + 1e-9)
    assert n_missing == int(outside.sum()) == 23 and n_total == 144
    assert n_clamped == 0
    assert vals.min() >= 1.0 and vals.max() <= 11.0
    assert np.allclose(vals[~outside], 1 + 10 * xy[~outside, 0], atol=1e-12)


def test_the_default_tolerance_extrapolates_and_the_transfer_does_not():
    r"""The raw cross-mesh interpolate, at Firedrake's default tolerance,
    continues 1 + x past the outline (0.97 at x = -0.03, below the source
    minimum of 1): the 2 km Budd snapshot onto the buffered 5 km mesh gave
    a prior in [-164.7, 921.1] from [1.0, 783.7] (issue 81). The transfer
    treats that band as outside the source."""
    source, target = _shifted_meshes()
    x, _y = fd.SpatialCoordinate(source)
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(1.0 + x)
    raw = fd.Function(fd.FunctionSpace(target, "CG", 1)).interpolate(
        f, allow_missing_dofs=True, default_missing_val=np.nan
    )
    assert np.nanmin(raw.dat.data_ro) < 1.0 - 1e-6
    g = fd.Function(fd.FunctionSpace(target, "CG", 1))
    n_missing, n_total, n_clamped = interpolate_with_fill(g, f, 1.0)
    xy = target.coordinates.dat.data_ro
    truly_outside = (xy < -1e-9).any(axis=1) | (xy > 1.0 + 1e-9).any(axis=1)
    vals = g.dat.data_ro
    assert n_clamped == 0
    assert n_missing == int(truly_outside.sum()) and n_total == 169
    assert np.all(vals[truly_outside] == 1.0)
    assert vals.min() >= 1.0 and vals.max() <= 2.0 + 1e-12
    assert np.allclose(vals[~truly_outside], 1.0 + xy[~truly_outside, 0], atol=1e-9)


def test_the_source_tolerance_is_restored_after_a_transfer():
    source, target = _meshes()
    source.tolerance = 0.25
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(fd.Constant(1.0))
    g = fd.Function(fd.FunctionSpace(target, "CG", 1))
    interpolate_with_fill(g, f, 7.0)
    assert source.tolerance == 0.25
    with pytest.raises(RuntimeError, match="inside"):
        with strict_location(source):
            assert source.tolerance == STRICT_TOLERANCE
            raise RuntimeError("inside")
    assert source.tolerance == 0.25
