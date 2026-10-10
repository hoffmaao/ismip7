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
        "ISMIP7_INVERSION_KSP_RTOL",
        "ISMIP7_KSP_RTOL",
        "ISMIP7_CONDENSED_KSP_ATOL_FACTOR",
        "ISMIP7_CONDENSED_PETSC_OPTIONS",
        "ISMIP7_DIRECT_FORWARD_MAXIT",
        "ISMIP7_DIRECT_FORWARD_DTOL",
        "ISMIP7_FINAL_SNES_STOL",
        "ISMIP7_FINAL_SNES_MAXIT",
    ):
        monkeypatch.delenv(name, raising=False)


def test_default_is_scpc_gamg_whatever_the_lane_solver(monkeypatch):
    assert sc.inversion_solver_mode() == "scpc_gamg"
    # The timing campaign exports scpc_mumps for the lane; the taped solve
    # has its own knob and stays put.
    monkeypatch.setenv("ISMIP7_DIAGNOSTIC_LINEAR_SOLVER", "scpc_mumps")
    assert sc.inversion_solver_mode() == "scpc_gamg"
    monkeypatch.setenv("ISMIP7_INVERSION_LINEAR_SOLVER", " FULL_MUMPS ")
    assert sc.inversion_solver_mode() == "full_mumps"


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
def test_scpc_modes_are_the_forward_options_at_a_tighter_krylov_tolerance(mode):
    params = sc.inversion_state_parameters(mode)
    forward = sc.diagnostic_solver_parameters(mode)
    assert params["mat_type"] == "matfree"
    assert params["pc_python_type"].endswith("ISMIP7SCPC")
    # the shared NLEQ-ERR, the transient's line search
    assert params["snes_linesearch_type"] == forward["snes_linesearch_type"] == "nleqerr"
    assert params["ksp_rtol"] == 1e-8
    assert forward["ksp_rtol"] == 1e-6
    if mode == "scpc_gamg":
        # the condensed solve's absolute tolerance follows the outer one
        assert params["condensed_field_ksp_atol"] == pytest.approx(0.5e-8)
        assert forward["condensed_field_ksp_atol"] == pytest.approx(0.5e-6)
    for key in ("ksp_rtol", "condensed_field_ksp_atol"):
        params.pop(key, None)
        forward.pop(key, None)
    assert params == forward


def test_a_rungs_condensed_options_still_come_last(monkeypatch):
    monkeypatch.setenv("ISMIP7_CONDENSED_PETSC_OPTIONS", "ksp_atol=1e-9 pc_gamg_threshold=0.02")
    params = sc.inversion_state_parameters("scpc_gamg")
    assert params["condensed_field_ksp_atol"] == "1e-9"
    assert params["condensed_field_pc_gamg_threshold"] == "0.02"


def test_the_inversion_knobs_override_and_leave_the_transient_alone(monkeypatch):
    monkeypatch.setenv("ISMIP7_INVERSION_SNES_LINESEARCH", "bt")
    monkeypatch.setenv("ISMIP7_INVERSION_KSP_RTOL", "1e-10")
    params = sc.inversion_state_parameters("scpc_gamg")
    assert params["snes_linesearch_type"] == "bt"
    assert params["ksp_rtol"] == 1e-10
    assert params["condensed_field_ksp_atol"] == pytest.approx(0.5e-10)
    forward = sc.diagnostic_solver_parameters("scpc_gamg")
    assert forward["snes_linesearch_type"] == "nleqerr"
    assert forward["ksp_rtol"] == 1e-6
    # full_mumps is the reference: neither knob reaches it
    assert sc.inversion_state_parameters("full_mumps")["snes_linesearch_type"] == "nleqerr"
    assert "ksp_rtol" not in sc.inversion_state_parameters("full_mumps")


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


@pytest.mark.parametrize("mode", sc.INVERSION_SOLVER_MODES)
def test_the_direct_solve_keeps_the_mode_and_bounds_the_trial(monkeypatch, mode):
    # PR 158's direct forward under every inversion solver: the mode's
    # linear options, a relative test only, a live step-size exit, and a
    # lost trial after 30 Newton iterations or a millionfold growth
    params = sc.inversion_state_parameters(mode)
    direct = sc.direct_forward_parameters(params)
    assert "snes_atol" in params and "snes_atol" not in direct
    assert direct["snes_max_it"] == 30
    assert direct["snes_divergence_tolerance"] == 1e6
    assert direct["snes_stol"] == float(sc.FINAL_SNES_STOL_DEFAULT)
    for key, value in params.items():
        if key not in ("snes_atol", "snes_max_it", "snes_stol",
                       "snes_divergence_tolerance"):
            assert direct[key] == value, key
    monkeypatch.setenv("ISMIP7_DIRECT_FORWARD_MAXIT", "12")
    monkeypatch.setenv("ISMIP7_DIRECT_FORWARD_DTOL", "1e4")
    direct = sc.direct_forward_parameters(params)
    assert (direct["snes_max_it"], direct["snes_divergence_tolerance"]) == (12, 1e4)
