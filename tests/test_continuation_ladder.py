r"""The exponent continuation restores its guess and re-ramps with more steps
when a rung diverges, and gives up only when every rung has."""
import numpy as np
import pytest

firedrake = pytest.importorskip("firedrake")
from firedrake import Constant, Function, FunctionSpace, UnitSquareMesh  # noqa: E402
from firedrake.exceptions import ConvergenceError  # noqa: E402

from icepack2_tools.continuation import ladder, ramp_exponents  # noqa: E402


def test_the_ladder_doubles_twice():
    assert ladder(8) == (8, 16, 32)
    assert ladder("5") == (5, 10, 20)
    with pytest.raises(ValueError):
        ladder(0)


@pytest.fixture
def state():
    mesh = UnitSquareMesh(2, 2)
    z = Function(FunctionSpace(mesh, "CG", 1)).assign(0.1)
    return z, Constant(1.0), Constant(1.0)


def _recording_solver(z, n_flow, m_slide, diverge_at=()):
    r"""A solve that advances ``z`` by the step's parameter and raises at the
    ``(attempt, step)`` pairs in ``diverge_at``, recording every call."""
    calls = []

    def solve(attempt, step, steps, t):
        calls.append((attempt, step, steps, round(t, 6),
                      round(float(n_flow), 6), round(float(m_slide), 6)))
        z.assign(z + Constant(t))
        if (attempt, step) in diverge_at:
            raise ConvergenceError("diverged")

    return solve, calls


def test_a_clean_ramp_takes_the_first_rung_to_the_targets(state):
    z, n_flow, m_slide = state
    solve, calls = _recording_solver(z, n_flow, m_slide)
    reports = []
    assert ramp_exponents(solve, z, n_flow, m_slide, 3.0, 3.0, ladder(4),
                          report=reports.append) == (1, 4)
    assert [c[:3] for c in calls] == [(1, 1, 4), (1, 2, 4), (1, 3, 4), (1, 4, 4)]
    # exponents ramp linearly and together, ending at the targets
    assert [c[4] for c in calls] == [1.0, 1.666667, 2.333333, 3.0]
    assert [c[5] for c in calls] == [1.0, 1.666667, 2.333333, 3.0]
    assert float(n_flow) == 3.0 and float(m_slide) == 3.0
    assert reports == []


def test_a_divergence_restores_the_guess_and_doubles_the_steps(state):
    z, n_flow, m_slide = state
    solve, calls = _recording_solver(z, n_flow, m_slide, diverge_at={(1, 3)})
    reports = []
    assert ramp_exponents(solve, z, n_flow, m_slide, 3.0, 2.0, ladder(4),
                          report=reports.append) == (2, 8)
    first, second = [c for c in calls if c[0] == 1], [c for c in calls if c[0] == 2]
    assert len(first) == 3 and len(second) == 8
    # the second rung starts over from the entry guess and from exponents 1
    assert second[0][3:] == (0.0, 1.0, 1.0)
    # z carries only the second rung's increments: sum of t over 8 steps = 4
    assert np.allclose(z.dat.data_ro, 0.1 + 4.0)
    assert float(n_flow) == 3.0 and float(m_slide) == 2.0
    assert reports == ["  Continuation diverged at 4 steps; restarting with 8..."]


def test_every_rung_diverging_raises_and_leaves_the_entry_state(state):
    z, n_flow, m_slide = state
    solve, calls = _recording_solver(
        z, n_flow, m_slide, diverge_at={(1, 1), (2, 1), (3, 1)})
    reports = []
    with pytest.raises(ConvergenceError):
        ramp_exponents(solve, z, n_flow, m_slide, 3.0, 3.0, ladder(2),
                       report=reports.append)
    assert [c[:3] for c in calls] == [(1, 1, 2), (2, 1, 4), (3, 1, 8)]
    assert np.allclose(z.dat.data_ro, 0.1)
    assert float(n_flow) == 1.0 and float(m_slide) == 1.0
    assert reports[-1] == "  Continuation diverged at 8 steps - giving up."
