r"""A 2 km MAP can warm-start an inversion on the 1 km mesh: the helpers the
inversion uses to decide what the warm start supplies."""
import numpy as np
import pytest

firedrake = pytest.importorskip("firedrake")
from firedrake import Function, FunctionSpace, UnitSquareMesh  # noqa: E402

from icepack2_tools.transfer import interpolate_with_fill, meshes_match  # noqa: E402


def test_the_same_mesh_matches_and_a_finer_one_does_not():
    coarse, fine = UnitSquareMesh(4, 4), UnitSquareMesh(8, 8)
    assert meshes_match(coarse, coarse)
    assert meshes_match(coarse, UnitSquareMesh(4, 4))
    assert not meshes_match(coarse, fine)


def test_a_warm_start_from_a_smaller_mesh_fills_the_prior_with_cold_ice():
    r"""Target dofs the source mesh does not cover take the stated fill (the
    inversion passes the forward's cold-ice baseline for the fluidity prior),
    so A = A_prior exp(phi) stays positive there."""
    source_mesh = UnitSquareMesh(4, 4)
    target_mesh = firedrake.RectangleMesh(8, 8, 1.5, 1.5)     # a wider domain
    src = Function(FunctionSpace(source_mesh, "CG", 1)).assign(40.0)
    tgt = Function(FunctionSpace(target_mesh, "CG", 1))
    n_missing, n_total, n_clamped = interpolate_with_fill(tgt, src, 1.0)
    assert 0 < n_missing < n_total and n_clamped == 0
    vals = tgt.dat.data_ro
    assert set(np.round(vals, 9)) == {1.0, 40.0}
    assert vals.min() > 0.0
