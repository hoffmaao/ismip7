r"""ISMIP7_INVERSION_LINEAR_SOLVER: the options of the inversion's taped
forward and of the adjoint tlm_adjoint solves against it. Pure option
dictionaries, no Firedrake."""

import pytest

from icepack2_tools import solverconfig as sc


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for name in (
        "ISMIP7_INVERSION_LINEAR_SOLVER",
        "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER",
        "ISMIP7_SNES_ATOL",
        "ISMIP7_SNES_RTOL",
        "ISMIP7_SNES_LINESEARCH",
        "ISMIP7_INVERSION_SNES_LINESEARCH",
    ):
        monkeypatch.delenv(name, raising=False)


def test_default_is_full_mumps_whatever_the_lane_solver(monkeypatch):
    assert sc.inversion_solver_mode() == "full_mumps"
    # The timing campaign exports scpc_mumps for the lane; the taped solve
    # has its own knob and stays put.
    monkeypatch.setenv("ISMIP7_DIAGNOSTIC_LINEAR_SOLVER", "scpc_mumps")
    assert sc.inversion_solver_mode() == "full_mumps"
    monkeypatch.setenv("ISMIP7_INVERSION_LINEAR_SOLVER", " SCPC_GAMG ")
    assert sc.inversion_solver_mode() == "scpc_gamg"


@pytest.mark.parametrize("mode", ["schur_gamg", "schur_mumps", "mumps", "lu"])
def test_unqualified_modes_are_refused(monkeypatch, mode):
    monkeypatch.setenv("ISMIP7_INVERSION_LINEAR_SOLVER", mode)
    with pytest.raises(ValueError, match="ISMIP7_INVERSION_LINEAR_SOLVER"):
        sc.inversion_solver_mode()


def test_full_mumps_is_the_inversion_reference_unchanged():
    # The dictionary the driver hardcoded before the knob existed; a MAP's
    # state_solver_parameters stamp is this, plus the forward's tolerance.
    expected = sc.nonlinear_solver_options()
    expected.update({
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 400,
        "mat_mumps_icntl_24": 1,
        "mat_mumps_cntl_3": 1e-12,
        "mat_mumps_icntl_4": 1,
    })
    assert sc.inversion_state_parameters("full_mumps") == expected


@pytest.mark.parametrize("mode", ["scpc_mumps", "scpc_gamg"])
def test_scpc_modes_are_the_forward_options(mode):
    params = sc.inversion_state_parameters(mode)
    forward = sc.diagnostic_solver_parameters(mode)
    assert params["mat_type"] == "matfree"
    assert params["pc_python_type"].endswith("ISMIP7SCPC")
    params.pop("snes_linesearch_type")
    forward.pop("snes_linesearch_type")
    assert params == forward


def test_gamg_backtracks_on_the_residual_and_exact_modes_do_not(monkeypatch):
    monkeypatch.delenv("ISMIP7_INVERSION_SNES_LINESEARCH", raising=False)
    monkeypatch.delenv("ISMIP7_SNES_LINESEARCH", raising=False)
    shared = sc.nonlinear_solver_options()["snes_linesearch_type"]
    assert shared == "nleqerr"
    assert sc.inversion_state_parameters("scpc_gamg")["snes_linesearch_type"] == "bt"
    assert sc.inversion_state_parameters("scpc_mumps")["snes_linesearch_type"] == shared
    assert sc.inversion_state_parameters("full_mumps")["snes_linesearch_type"] == shared
    # the transient's scpc_gamg is untouched
    assert sc.diagnostic_solver_parameters("scpc_gamg")["snes_linesearch_type"] == shared
    monkeypatch.setenv("ISMIP7_INVERSION_SNES_LINESEARCH", "nleqerr")
    assert sc.inversion_state_parameters("scpc_gamg")["snes_linesearch_type"] == "nleqerr"


@pytest.mark.parametrize("mode", ["full_mumps", "scpc_mumps", "scpc_gamg"])
def test_adjoint_is_one_linear_solve_with_no_absolute_exit(mode):
    params = sc.inversion_state_parameters(mode)
    params["snes_atol"] = 1e-3
    params["snes_monitor"] = None
    adjoint = sc.inversion_adjoint_parameters(params)
    snes = {k: v for k, v in adjoint.items() if k.startswith("snes_")}
    if params.get("mat_type") == "matfree":
        # tlm_adjoint hands these to a LinearVariationalSolver, whose ksponly
        # default the forward's newtonls would otherwise replace
        assert snes == {"snes_type": "ksponly"}
    else:
        assert snes == {}
    # job 10432790: an absolute exit returned a zero adjoint
    assert "ksp_atol" not in adjoint
    # everything else is the forward's
    rest = {k: v for k, v in params.items() if not k.startswith("snes_")}
    assert {k: v for k, v in adjoint.items() if not k.startswith("snes_")} == rest
    # and the forward's dictionary is left alone
    assert params["snes_atol"] == 1e-3
