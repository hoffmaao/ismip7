r"""A MAP transferred onto a mesh with an ocean buffer must not be
extrapolated into the buffer (22 September 2026: fluidity prior
[-164.7, 921.1] from a source range [1.0, 783.7], then a singular local
block in the condensed solve). The transfer locates strictly, fills outside
the source and stays within the source range."""
import numpy as np
from firedrake import (Function, FunctionSpace, RectangleMesh, as_vector,
                       SpatialCoordinate, UnitSquareMesh, VectorFunctionSpace)

from icepack2_tools.mpi_stats import global_count
from icepack2_tools.transfer import outside_source, strict_transfer


def _meshes(shift=0.03):
    source = UnitSquareMesh(8, 8)
    # The target reaches past the source on every side, and its first row of
    # dofs sits ``shift`` outside the outline: closer than half a source cell
    # (0.0625 in reference L1 distance), which Firedrake's default tolerance
    # accepts as inside, and further than the strict tolerance rejects.
    target = RectangleMesh(12, 12, 1.5, 1.5)
    target.coordinates.dat.data[:] -= shift
    return source, target


def _source_field(source):
    x, y = SpatialCoordinate(source)
    return Function(FunctionSpace(source, "CG", 1), name="f").interpolate(1.0 + x)


def test_the_default_tolerance_extrapolates():
    source, target = _meshes()
    f = _source_field(source)
    g = Function(FunctionSpace(target, "CG", 1)).interpolate(
        f, allow_missing_dofs=True, default_missing_val=0.0)
    # Firedrake's 0.5 reference-cell tolerance lets dofs just outside the
    # source be located and the linear field continued past its range:
    # 1 + x at x = -0.03 reads 0.97, below the source minimum of 1.
    assert g.dat.data_ro.min() < 1.0 - 1e-6


def test_strict_transfer_fills_and_bounds():
    source, target = _meshes()
    f = _source_field(source)
    g, info = strict_transfer(f, FunctionSpace(target, "CG", 1), fill=1.0)
    data = g.dat.data_ro
    assert 1.0 - 1e-12 <= data.min() and data.max() <= 2.0 + 1e-12
    outside = outside_source(source, g.function_space())
    xy = Function(VectorFunctionSpace(target, "CG", 1)).interpolate(
        SpatialCoordinate(target)).dat.data_ro
    truly_outside = (xy < -1e-9).any(axis=1) | (xy > 1.0 + 1e-9).any(axis=1)
    assert np.array_equal(outside, truly_outside)
    assert np.all(data[outside] == 1.0)
    assert np.allclose(data[~outside], 1.0 + xy[~outside, 0], atol=1e-9)
    # n_outside is a global dof count; the mask here is this rank's slice.
    assert info["n_outside"] == global_count(truly_outside, target.comm)
    assert info["source_range"] == (1.0, 2.0)


def test_a_vector_field_is_filled_row_wise():
    source, target = _meshes()
    x, y = SpatialCoordinate(source)
    u = Function(VectorFunctionSpace(source, "CG", 1), name="u").interpolate(
        as_vector((x, -y)))
    g, info = strict_transfer(u, VectorFunctionSpace(target, "CG", 1), fill=0.0)
    outside = outside_source(source, FunctionSpace(target, "CG", 1))
    assert np.all(g.dat.data_ro[outside] == 0.0)
    assert info["n_outside"] > 0


def test_the_mask_is_strict_before_any_transfer():
    source, target = _meshes()
    space = FunctionSpace(target, "CG", 1)
    outside = outside_source(source, space)
    xy = Function(VectorFunctionSpace(target, "CG", 1)).interpolate(
        SpatialCoordinate(target)).dat.data_ro
    truly_outside = (xy < -1e-9).any(axis=1) | (xy > 1.0 + 1e-9).any(axis=1)
    assert np.array_equal(outside, truly_outside)


def test_a_fill_outside_the_source_range_is_reported_clipped():
    source, target = _meshes()
    f = _source_field(source)
    g, info = strict_transfer(f, FunctionSpace(target, "CG", 1), fill=0.0)
    outside = outside_source(source, g.function_space())
    assert info["fill"] == 1.0
    assert np.all(g.dat.data_ro[outside] == 1.0)
