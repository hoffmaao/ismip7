"""MUMPS workspace margin: one definition, under the name MUMPS reads.

Firedrake's own default margin is written under ``pc_factor_`` and never
reaches MUMPS (firedrake#5311), so every project solver that factors with
MUMPS must name ``mat_mumps_icntl_14`` itself. At MUMPS's default the
Greenland 1 km mass-matrix Cholesky failed 13 of 150 factorizations on 32
ranks (INFOG(1) = -9) and ended gris2 inversion links."""
import pytest

from icepack2_tools import solverconfig as sc

FACTOR_PCS = {"lu", "cholesky"}


def test_margin_is_named_as_mumps_reads_it():
    assert sc.mumps_workspace() == {"mat_mumps_icntl_14": sc.MUMPS_WORKSPACE_PERCENT}
    assert sc.mumps_workspace("condensed_field_") == {
        "condensed_field_mat_mumps_icntl_14": sc.MUMPS_WORKSPACE_PERCENT}
    for pc_type in FACTOR_PCS:
        params = sc.direct_solver_parameters(pc_type)
        assert params["pc_type"] == pc_type
        assert params["mat_mumps_icntl_14"] == sc.MUMPS_WORKSPACE_PERCENT
        assert "ksp_error_if_not_converged" in params
        assert not any(key.startswith("pc_factor_mat_mumps") for key in params)
    with pytest.raises(ValueError):
        sc.direct_solver_parameters("ilu")


@pytest.mark.parametrize("mode,prefix", [
    ("scpc_mumps", "condensed_field_"), ("schur_mumps", "fieldsplit_1_")])
def test_condensed_factorizations_carry_the_margin(monkeypatch, mode, prefix):
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: False)
    params = sc.diagnostic_solver_parameters(mode)
    assert params[f"{prefix}mat_mumps_icntl_14"] == sc.MUMPS_WORKSPACE_PERCENT


def test_full_jacobian_keeps_its_margin(monkeypatch):
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: False)
    full = sc.MUMPS_FULL_JACOBIAN_WORKSPACE_PERCENT
    assert sc.diagnostic_solver_parameters("full_mumps")["mat_mumps_icntl_14"] == full
    assert sc.inversion_state_parameters("full_mumps")["mat_mumps_icntl_14"] == full


def test_prior_uses_the_shared_definition():
    from icepack2_tools import prior
    assert prior.MASS_CHOLESKY == sc.direct_solver_parameters("cholesky")


# --- live: the margin reaches the factor matrix ------------------------------

fd = pytest.importorskip("firedrake")


def _mass(n=8):
    mesh = fd.UnitSquareMesh(n, n)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    u, v = fd.TrialFunction(Q), fd.TestFunction(Q)
    return Q, fd.assemble(fd.inner(u, v) * fd.dx), fd.assemble(fd.inner(fd.Constant(1.0), v) * fd.dx)


@pytest.mark.parametrize("pc_type", sorted(FACTOR_PCS))
def test_margin_reaches_the_factor_matrix(pc_type):
    Q, M, b = _mass()
    solver = fd.LinearSolver(M, solver_parameters=sc.direct_solver_parameters(pc_type))
    x = fd.Function(Q)
    solver.solve(x, b)
    factor = solver.ksp.getPC().getFactorMatrix()
    assert factor.getMumpsIcntl(14) == sc.MUMPS_WORKSPACE_PERCENT
    assert abs(x.dat.data_ro - 1.0).max() < 1e-10


def test_firedrake_default_margin_is_lost():
    """The reason for all of this (firedrake#5311): with Firedrake's defaults
    MUMPS runs at its own margin. If this fails, Firedrake fixed it."""
    Q, M, b = _mass()
    solver = fd.LinearSolver(M)
    solver.solve(fd.Function(Q), b)
    factor = solver.ksp.getPC().getFactorMatrix()
    assert factor.getMumpsIcntl(14) != 200


def test_failed_factorization_raises_at_its_source():
    from firedrake.petsc import PETSc
    mesh = fd.UnitSquareMesh(4, 4)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    u, v = fd.TrialFunction(Q), fd.TestFunction(Q)
    Z = fd.assemble(fd.Constant(0.0) * fd.inner(u, v) * fd.dx)   # singular
    b = fd.assemble(fd.inner(fd.Constant(1.0), v) * fd.dx)
    solver = fd.LinearSolver(Z, solver_parameters=sc.direct_solver_parameters("lu"))
    with pytest.raises(PETSc.Error):
        solver.solve(fd.Function(Q), b)
