# Uncertainty from the inferred friction and fluidity

**Goal.** The inversion returns one state for the basal friction and the
ice-shelf fluidity (the rheology field that absorbs fabric, damage and
temperature). This work asks how much of that state the observations
actually pin down, and how much the rest moves the dynamics the projections
start from. The first step, here, is the Laplace posterior at the MAP and the
spread it puts on the ice discharged from grounded ice at the start of the
projections, per drainage basin.

## Method

At the MAP the posterior is approximated as a Gaussian whose covariance is
the prior covariance `Gamma` updated by the Gauss-Newton Hessian `H` of the
data term, `Gamma_post = (H + Gamma^-1)^-1`, following Recinos et al. (2023).
`H = J^T R^-1 J`, with `J` the derivative of the momentum balance's velocity
with respect to the controls and `R^-1` the velocity term's own weights.
`antarctica/scripts/run_eigendec.py` computes the leading generalised
eigenpairs `H v_i = lambda_i Gamma^-1 v_i`. A mode with `lambda_i >> 1` is
set by the data; one with `lambda_i << 1` is known no better than the prior.
In that basis

    Gamma_post = Gamma - sum_i lambda_i / (1 + lambda_i) v_i v_i^T.

A quantity `Q` the model computes from the controls, linearised at the MAP
with gradient `g`, then has posterior variance

    sigma_post^2 = g . Gamma g - sum_i lambda_i / (1 + lambda_i) (g . v_i)^2

(`icepack2_tools.laplace.posterior_variance`). Leaving out modes leaves out
only reductions, so a truncated set of modes bounds the variance from above.

The quantities are the ice leaving grounded cells at the MAP state: the
upwind flux across the grounding line into floating cells, signed as the
ISMIP7 writer books `ligroundf`, plus the outflow of grounded cells across
the mesh boundary (`icepack2_tools.laplace.grounded_discharge_forms`). As in
the writer, every cell on a bed at or above sea level counts as grounded
(`icepack2_tools.laplace.grounded_indicator`), so ice flowing onto bare land
books nowhere. The totals are in
Gt/yr, in total, per IMBIE2 basin, and the grounding-line flux alone. Each
needs one adjoint solve of the momentum balance on the same recorded solve
the Hessian uses.

## What is and is not included

- Controls: the friction control `alpha` (square-root form) and the log
  fluidity `phi` on floating ice, with the bi-Laplacian priors the MAP
  records.
- Likelihood: the velocity term only. The production MAP's objective also
  scored dH/dt (`dhdt_weight` 1), so the posterior here is wider than the
  inversion's: the spreads are upper bounds. The output records the MAP's
  `dhdt_weight` beside them.
- Geometry: fixed (BedMachine after the relaxation). Thickness is not a
  control of the inversion, so its uncertainty is not part of this posterior.
- Time: the quantities are snapshots of the momentum balance at the MAP, the
  state the projections start from. Their spread is a first-order measure of
  how uncertain the initial dynamics are; how it grows through the projection
  time loop needs one of the next steps below.

## Running it

    ISMIP7_INVERSION=<map.h5> mpiexec -n 32 python antarctica/scripts/run_eigendec.py \
        --modes 40 --maxiter 2 --tag <tag>

under the environment the MAP was inverted in. The gradients and prior
spreads are printed before the eigensolve; the eigendecomposition file,
`eigenvalues_<lc>_<tag>.txt` and `uncertainty_<lc>_<tag>.json` (value, prior
and posterior sigma per quantity, and each mode's share of the reduction) are
written beside the results. Run once with `--check` first: it compares the
total discharge's adjoint gradient with a central difference of two untaped
solves, and the Gauss-Newton action along that direction with the norm of the
velocity change two more solves give, times one Gauss-Newton action, then
stops (the checks move the controls, so they cannot share a run with the
eigensolve).

Size the eigensolve to the job's wall-clock limit from that action cost.
ARPACK takes one Hessian action per Lanczos step, at most
`1 + ncv + maxiter (ncv - modes)` with `ncv = 2 modes + 1` (`--ncv` changes
it), so `--maxiter` bounds the run: `--modes 40 --maxiter 2` is at most 164
actions, nine hours at the 190 s per action the production MAP takes on 32
ranks (2 km; the 40-mode run needed 103 of them). When the restarts run out
first, the modes ARPACK has converged are kept and the json records
`converged: false`; the spreads they give are still upper bounds. `--tol`
(1e-6) is the relative accuracy asked of the eigenvalues, far finer than the
truncation.

Memory sets the other limit. At 2 km the eigensolve holds about 2.3 GB per
rank between actions and peaks near 1.4 times that during one, so 32 ranks
fit a 240 GB node; the log prints the resident set after every action, and
it must stay level. It did not before `patch_tangent_linear_cleanup`
(`icepack2_tools/laplace.py`): tlm_adjoint 0.9.1dev leaves each action's
tangent-linear copies of the forward variables attached to them, 0.6 GB per
action here, which ended a 100-mode run out of memory after 148 actions.

## Next steps

1. The dH/dt term in the Hessian, through the inversion's own objective, so
   the posterior uses every datum the MAP did.
2. Samples of the posterior run through the forward suite, for the spread
   of the projections themselves under the production physics.
