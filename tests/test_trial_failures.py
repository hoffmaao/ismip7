r"""``optimization.TrialFailures`` against scipy's L-BFGS-B, which the
inversion's scipy path runs with ``gtol=0``.

A failed forward or adjoint at a line-search trial returns ten times the last
good objective and a zero gradient. Issue 153's RC link 11461566 hit three in
a row at evaluation 107 and L-BFGS-B stopped with a CONVERGENCE message at a
gradient norm of 1.03; the chain then marked the MAP done. The tracker has to
tell that stop apart from a minimisation that met failed trials along the way
and recovered.
"""
import numpy as np
from scipy.optimize import minimize

from icepack2_tools.optimization import TrialFailures

HESSIAN = np.diag(np.linspace(1.0, 50.0, 20))


def _minimise(fails, regularisation=0.0):
    r"""L-BFGS-B on a quadratic, as the driver calls it; ``fails(n, x, x_good)``
    says whether evaluation ``n`` at ``x`` fails. ``regularisation`` is a
    constant added to the objective, so the total the failed trial scales from
    can dwarf the quadratic."""
    tracker = TrialFailures()
    good = {"f": None, "x": None, "n": 0}

    def objective_and_gradient(x):
        good["n"] += 1
        if good["x"] is not None and fails(good["n"], x, good["x"]):
            return tracker.failed(good["f"], x.size)
        f = 0.5 * x @ HESSIAN @ x + regularisation
        good.update(f=f, x=x.copy())
        tracker.succeeded()
        return f, HESSIAN @ x

    result = minimize(objective_and_gradient, np.ones(20), jac=True,
                      method="L-BFGS-B",
                      options={"maxiter": 300, "ftol": 1e-12, "gtol": 0})
    return result, tracker, float(np.linalg.norm(HESSIAN @ result.x)), good["x"]


def test_a_stop_after_failed_trials_is_flagged():
    r"""Every trial fails from evaluation 13 on, as when the direct forward
    stalls at the residual floor near the accepted point. scipy calls the
    stop a success; the gradient there is far from small."""
    result, tracker, grad_norm, _x_good = _minimise(lambda n, x, x_good: n > 12)
    assert result.success
    assert grad_norm > 0.1
    assert tracker.stopped_on_failures
    assert tracker.trailing == tracker.total > 0


def test_failed_trials_the_line_search_recovers_from_are_not_flagged():
    r"""Trials that step too far fail and the line search backtracks past
    them; the minimisation converges and the tracker leaves it final."""
    result, tracker, grad_norm, _x_good = _minimise(
        lambda n, x, x_good: np.linalg.norm(x - x_good) > 0.05)
    assert tracker.total > 0
    assert grad_norm < 1e-4
    assert not tracker.stopped_on_failures


def test_a_minimisation_without_failures_is_not_flagged():
    result, tracker, grad_norm, _x_good = _minimise(lambda n, x, x_good: False)
    assert grad_norm < 1e-4
    assert tracker.total == 0
    assert not tracker.stopped_on_failures


def test_a_stop_after_failed_trials_holds_the_last_good_point_under_heavy_regularisation():
    r"""Regularisation a hundred times the misfit at the start: scaled from the
    total, every failed trial stays above the current objective, L-BFGS-B
    rejects it and the result is the last point whose forward succeeded."""
    result, tracker, _grad_norm, x_good = _minimise(
        lambda n, x, x_good: n > 12, regularisation=100 * 0.5 * np.trace(HESSIAN))
    assert tracker.stopped_on_failures
    assert np.array_equal(result.x, x_good)


def test_a_failed_trial_returns_ten_times_the_last_objective_and_no_gradient():
    objective, gradient = TrialFailures().failed(3.0, 4)
    assert objective == 30.0
    assert np.array_equal(gradient, np.zeros(4))
