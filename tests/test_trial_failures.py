r"""``optimization.TrialFailures`` against scipy's L-BFGS-B, which the
inversion's scipy path runs with ``gtol=0``.

A failed forward or adjoint at a line-search trial returns ten times the last
good objective and a zero gradient. Issue 153's RC link 11461566 hit three in
a row at evaluation 107 and L-BFGS-B stopped with a CONVERGENCE message at a
gradient norm of 1.03; the chain then marked the MAP done. The tracker has to
tell that stop apart from a minimisation that met failed trials along the way
and recovered.

``optimization.AcceptedIterate`` keeps the evaluation scipy returns, which is
the accepted iterate and rarely the last evaluation: the publishing solve and
the MAP's record need that evaluation's state and objective.
"""
import numpy as np
from scipy.optimize import minimize

from icepack2_tools.optimization import AcceptedIterate, TrialFailures

HESSIAN = np.diag(np.linspace(1.0, 50.0, 20))


def _minimise(fails, regularisation=0.0, scale_from="total"):
    r"""L-BFGS-B on a quadratic misfit plus a constant ``regularisation``, as
    the driver calls it; ``fails(n, x, x_good)`` says whether evaluation ``n``
    at ``x`` fails, and a failed trial scales from the last good
    ``scale_from`` ("misfit" or "total")."""
    tracker = TrialFailures()
    good = {"misfit": None, "total": None, "x": None, "n": 0}

    def objective_and_gradient(x):
        good["n"] += 1
        if good["x"] is not None and fails(good["n"], x, good["x"]):
            return tracker.failed(good[scale_from], x.size)
        misfit = 0.5 * x @ HESSIAN @ x
        good.update(misfit=misfit, total=misfit + regularisation, x=x.copy())
        tracker.succeeded()
        return good["total"], HESSIAN @ x

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
    r"""Regularisation a hundred times the starting misfit and every trial
    failing from evaluation 13 on. Scaled from the misfit alone, a failed
    trial sits below the current total, L-BFGS-B accepts it and the result is
    a point whose forward failed; scaled from the total, every failed trial is
    rejected and the result is the last point whose forward succeeded."""
    regularisation = 100 * 0.5 * np.trace(HESSIAN)
    for scale_from, holds in (("misfit", False), ("total", True)):
        result, tracker, _grad_norm, x_good = _minimise(
            lambda n, x, x_good: n > 12, regularisation=regularisation,
            scale_from=scale_from)
        assert tracker.stopped_on_failures
        assert np.array_equal(result.x, x_good) == holds, scale_from


def test_a_failed_trial_returns_ten_times_the_last_objective_and_no_gradient():
    objective, gradient = TrialFailures().failed(3.0, 4)
    assert objective == 30.0
    assert np.array_equal(gradient, np.zeros(4))


def _minimise_against_a_wall(wall):
    r"""The quadratic misfit with every trial past ``x.sum() < wall`` failing,
    as a forward that stops converging past some controls; the minimiser at
    0 lies past the wall. Each good evaluation records its objective in an
    ``AcceptedIterate``, as the driver records its mixed state."""
    tracker, accepted = TrialFailures(), AcceptedIterate()
    good = {"total": None, "x": None}

    def objective_and_gradient(x):
        if good["x"] is not None and x.sum() < wall:
            return tracker.failed(good["total"], x.size)
        good.update(total=0.5 * x @ HESSIAN @ x, x=x.copy())
        accepted.evaluated(good["x"], good["total"])
        tracker.succeeded()
        return good["total"], HESSIAN @ x

    result = minimize(objective_and_gradient, np.ones(20), jac=True,
                      method="L-BFGS-B", callback=accepted.iterate,
                      options={"maxiter": 300, "ftol": 1e-12, "gtol": 0})
    return result, tracker, accepted, good["x"]


def test_the_returned_controls_are_the_accepted_evaluation_and_not_the_last():
    r"""Against the wall at 5 L-BFGS-B stops right after a failed trial, at
    2 it stops ABNORMAL with none trailing. Both times result.x is an earlier
    evaluation than the last good one, whose state the publishing solve used
    to start from; after the failure stop result.fun is the failed trial's
    inflated value. The record holds the accepted evaluation's objective."""
    for wall, stopped in ((5.0, True), (2.0, False)):
        result, tracker, accepted, x_last = _minimise_against_a_wall(wall)
        assert tracker.stopped_on_failures == stopped, wall
        assert not np.array_equal(result.x, x_last), wall
        objective = 0.5 * result.x @ HESSIAN @ result.x
        assert accepted.at(result.x) == objective, wall
        if stopped:
            assert result.fun > objective


def test_the_starting_point_is_accepted_until_a_callback_moves_it():
    accepted = AcceptedIterate()
    x0, x1 = np.zeros(3), np.ones(3)
    accepted.evaluated(x0, "start")
    accepted.evaluated(x1, "trial")
    assert accepted.at(x0) == "start"
    accepted.iterate(np.full(3, 2.0))
    assert accepted.at(x0) == "start"
    accepted.iterate(x1.copy())
    assert accepted.at(x1) == "trial"
    assert accepted.at(x0) is None
