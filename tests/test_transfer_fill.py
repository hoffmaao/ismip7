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


# ── Smooth extension (ISMIP7_TRANSFER_FILL=extend) ─────────────────────────

def _ring_meshes(n=32):
    # Target vertices on x <= 1 coincide with the source's, so the ring
    # 1 < x <= 1.5 is exactly the missing set.
    source = fd.UnitSquareMesh(n, n)
    target = fd.RectangleMesh(3 * n // 2, n, 1.5, 1.0)
    return source, target


def test_the_extension_is_the_harmonic_continuation_of_the_outline_values():
    r"""cos(pi y) on the outline x = 1 continues into the ring 1 < x < 1.5
    with zero normal derivative on its other sides as
    cos(pi y) cosh(pi (x - 1.5)) / cosh(pi / 2), which is harmonic."""
    source, target = _ring_meshes()
    _x, y = fd.SpatialCoordinate(source)
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(fd.cos(np.pi * y))
    g = fd.Function(fd.FunctionSpace(target, "CG", 1), name="theta")
    n_missing, _n_total, n_clamped = interpolate_with_fill(g, f, 0.0, extend=True)
    xy = target.coordinates.dat.data_ro
    ring = xy[:, 0] > 1.0 + 1e-9
    assert n_missing == int(ring.sum()) and n_clamped == 0
    exact = (np.cos(np.pi * xy[:, 1]) * np.cosh(np.pi * (xy[:, 0] - 1.5))
             / np.cosh(np.pi / 2))
    vals = g.dat.data_ro
    assert np.abs(vals[ring] - exact[ring]).max() < 2e-3
    assert np.allclose(vals[~ring], np.cos(np.pi * xy[~ring, 1]), atol=1e-12)


def _curvature_energy(g):
    # (K u)' M_lumped^-1 (K u): the bi-Laplacian's curvature term, with the
    # mass lumped so no solve is needed.
    V = g.function_space()
    u, v = fd.TrialFunction(V), fd.TestFunction(V)
    K = fd.assemble(fd.inner(fd.grad(u), fd.grad(v)) * fd.dx).petscmat
    m = fd.assemble(v * fd.dx)
    with g.dat.vec_ro as gv:
        ku = K.createVecLeft()
        K.mult(gv, ku)
        return float(np.sum(ku.array ** 2 / m.dat.data_ro))


def test_the_extension_removes_the_step_a_constant_fill_puts_at_the_outline():
    r"""The issue 153 case: controls of order one at the outline and zero in
    the ring. Measured against the source field's own energy, the extension
    adds 0.34 % of what the constant fill's step adds at n = 32, 0.08 % at
    n = 64: the step's cost grows as the mesh refines and the extension's
    does not."""
    source, target = _ring_meshes()
    x, y = fd.SpatialCoordinate(source)
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(
        2.0 + fd.sin(2 * np.pi * y) * x)
    Q = fd.FunctionSpace(target, "CG", 1)
    const, ext = fd.Function(Q), fd.Function(Q)
    interpolate_with_fill(const, f, 0.0)
    interpolate_with_fill(ext, f, 0.0, extend=True)
    e_const, e_ext, e_src = (_curvature_energy(g) for g in (const, ext, f))
    assert (e_ext - e_src) < 1e-2 * (e_const - e_src)


def test_a_log_extension_stays_positive_and_inside_the_source_range():
    source, target = _ring_meshes(16)
    _x, y = fd.SpatialCoordinate(source)
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(1.0 + 782.0 * y ** 4)
    g = fd.Function(fd.FunctionSpace(target, "CG", 1), name="fluidity_prior")
    interpolate_with_fill(g, f, 11.1, extend=True, log=True)
    xy = target.coordinates.dat.data_ro
    ring = xy[:, 0] > 1.0 + 1e-9
    vals = g.dat.data_ro
    assert vals.min() >= 1.0 and vals.max() <= 783.0 + 1e-9
    # the extension of log(A), not of A: at the far edge of the ring a
    # geometric mean of the outline values, well below their arithmetic mean
    assert np.all(vals[ring] > 0.0)
    far = ring & np.isclose(xy[:, 0], 1.5)
    assert np.median(vals[far]) < 0.5 * np.mean(1.0 + 782.0 * xy[far, 1] ** 4)


def test_a_region_the_source_never_reaches_keeps_the_fill():
    _source, target = _ring_meshes(8)
    Q = fd.FunctionSpace(target, "CG", 1)
    g = fd.Function(Q).assign(5.0)
    from icepack2_tools.transfer import harmonic_extension
    harmonic_extension(g, np.ones(g.dat.data_ro.shape[0], dtype=bool), 3.0)
    # the pure-Neumann block is pinned by a 1e-8 mass term, so rounding
    # in K x = 0 shows at 1e-8 relative
    assert np.allclose(g.dat.data_ro, 3.0, rtol=1e-6, atol=0)


def test_the_extension_refuses_what_it_cannot_extend():
    source, target = _ring_meshes(8)
    f = fd.Function(fd.FunctionSpace(source, "CG", 1)).interpolate(fd.Constant(1.0))
    g = fd.Function(fd.FunctionSpace(target, "CG", 1))
    with pytest.raises(ValueError, match="float fill"):
        interpolate_with_fill(g, f, fd.Function(g.function_space()), extend=True)
    f_dg = fd.Function(fd.FunctionSpace(source, "DG", 0)).interpolate(fd.Constant(1.0))
    g_dg = fd.Function(fd.FunctionSpace(target, "DG", 0))
    with pytest.raises(ValueError, match="scalar continuous"):
        interpolate_with_fill(g_dg, f_dg, 0.0, extend=True)
