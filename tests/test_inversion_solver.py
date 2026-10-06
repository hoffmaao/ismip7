"""The inversion's linear solver: static condensation by default, its adjoint options, the
PT-Scotch fallback, and the quadrature degree stamped into the residual."""
import sys
from pathlib import Path

import pytest

from icepack2_tools import solverconfig as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "antarctica" / "scripts"))


def test_condensation_is_the_inversion_default(monkeypatch):
    monkeypatch.delenv("ISMIP7_INVERSION_LINEAR_SOLVER", raising=False)
    assert sc.inversion_solver_mode() == "scpc_mumps"
    monkeypatch.setenv("ISMIP7_INVERSION_LINEAR_SOLVER", "full_mumps")
    assert sc.inversion_solver_mode() == "full_mumps"
    monkeypatch.setenv("ISMIP7_INVERSION_LINEAR_SOLVER", "lu_everything")
    with pytest.raises(ValueError):
        sc.inversion_solver_mode()


def test_adjoint_options_are_linear_matfree_and_tight(monkeypatch):
    monkeypatch.delenv("ISMIP7_ADJOINT_KSP_RTOL", raising=False)
    p = sc.adjoint_solver_parameters("scpc_mumps")
    assert not any(k.startswith("snes_") for k in p)
    assert p["mat_type"] == "matfree"        # reaches tlm_adjoint's matrix-free adjoint branch
    assert p["pc_python_type"] == "icepack2_tools.preconditioners.ISMIP7SCPC"
    assert p["ksp_rtol"] == 1e-10 and p["ksp_atol"] == 0.0
    monkeypatch.setenv("ISMIP7_ADJOINT_KSP_RTOL", "1e-12")
    assert sc.adjoint_solver_parameters("scpc_mumps")["ksp_rtol"] == 1e-12


def test_mumps_ordering_follows_the_petsc_build(monkeypatch):
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: False)
    opts = sc._mumps_options("x_")
    assert "x_mat_mumps_icntl_29" not in opts and opts["x_pc_factor_mat_solver_type"] == "mumps"
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: True)
    monkeypatch.delenv("ISMIP7_MUMPS_ANALYSIS", raising=False)
    opts = sc._mumps_options("x_")
    assert opts["x_mat_mumps_icntl_28"] == 2 and opts["x_mat_mumps_icntl_29"] == 1
    monkeypatch.setenv("ISMIP7_MUMPS_ANALYSIS", "sequential")
    assert "x_mat_mumps_icntl_28" not in sc._mumps_options("x_")


def test_quadrature_degree_reaches_every_derived_form():
    """The adjoint of the Jacobian of a stamped residual still carries the degree,
    so a solve that passes no form-compiler parameters integrates it the same way."""
    fd = pytest.importorskip("firedrake")
    import ufl
    from icepack2_tools.forms import with_quadrature_degree
    mesh = fd.UnitSquareMesh(2, 2)
    V = fd.FunctionSpace(mesh, "CG", 1)
    u, v = fd.Function(V), fd.TestFunction(V)
    F = with_quadrature_degree(u ** 3 * v * fd.dx + fd.avg(u) * fd.jump(v) * fd.dS, 4)
    for form in (F, fd.derivative(F, u), fd.adjoint(fd.derivative(F, u))):
        assert {itg.metadata()["quadrature_degree"] for itg in form.integrals()} == {4}
    kept = with_quadrature_degree(u * v * fd.dx(degree=7), 4)
    assert kept.integrals()[0].metadata()["quadrature_degree"] == 7
    assert isinstance(F, ufl.Form)


def test_mumps_analysis_leaves_the_cache_fingerprint_alone(monkeypatch):
    """The analysis follows the build and the knob; the record says which ran,
    and a cache prepared on one build is accepted on another."""
    import timing_campaign as tc
    monkeypatch.delenv("ISMIP7_MUMPS_ANALYSIS", raising=False)
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: True)
    parallel = sc.solver_provenance("scpc_mumps")
    assert parallel["diagnostic_petsc_options"]["condensed_field_mat_mumps_icntl_28"] == 2
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: False)
    no_ptscotch = sc.solver_provenance("scpc_mumps")
    assert "condensed_field_mat_mumps_icntl_28" not in no_ptscotch["diagnostic_petsc_options"]
    monkeypatch.setattr(sc, "_have_ptscotch", lambda: True)
    monkeypatch.setenv("ISMIP7_MUMPS_ANALYSIS", "sequential")
    sequential = sc.solver_provenance("scpc_mumps")
    assert "condensed_field_mat_mumps_icntl_28" not in sequential["diagnostic_petsc_options"]
    fingerprints = {tc.solver_configuration_fingerprint(p)
                    for p in (parallel, no_ptscotch, sequential)}
    assert len(fingerprints) == 1
    assert fingerprints != {tc.solver_configuration_fingerprint(
        sc.solver_provenance("full_mumps"))}


def test_mumps_analysis_knob_is_checked(monkeypatch):
    monkeypatch.setenv("ISMIP7_MUMPS_ANALYSIS", "distributed")
    with pytest.raises(ValueError):
        sc.mumps_analysis()
