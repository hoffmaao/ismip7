r"""Cross-mesh interpolation with a stated fill for the dofs the source mesh
does not cover (icepack2_tools.transfer.interpolate_with_fill).

The production case: a MAP inverted on a buffer-0 mesh transferred onto the
buffered production mesh. Every target dof in the 20 km ocean ring is
outside the source, and the fluidity prior must not become zero there.
"""
import numpy as np
import pytest

fd = pytest.importorskip("firedrake")

from icepack2_tools.transfer import interpolate_with_fill  # noqa: E402


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
    n_missing, n_total = interpolate_with_fill(g, f, 7.0)
    outside = _outside(target)
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
    n_missing, n_total = interpolate_with_fill(w, u, raster)
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
    n_missing, n_total = interpolate_with_fill(g, f, 7.0)
    assert n_missing == 0 and n_total == 25
    assert np.allclose(g.dat.data_ro, f.dat.data_ro)
