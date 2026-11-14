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
`H = J^T R^-1 J`, with `J` the derivative of what the data term observes
(the momentum balance's velocity, and under the MAP's own objective the
one-year thickness step and its filtered tendency) with respect to the
controls, and `R^-1` the data term's own weights.
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
- Likelihood: the MAP's own objective (`--likelihood map`, the default),
  rebuilt from the settings the MAP records: the velocity term with its
  dH/dt re-weighting (`dhdt_reweight`: up to 25 times the weight where the
  observed thinning reaches 1 m/yr on grounded ice), the dH/dt term through
  one upwind step of the MAP's transport and its Helmholtz filter
  (`dhdt_weight`, `dhdt_filter_l`, sigma 0.1 m/yr from the environment, as in
  the inversion), and the net term (`dhdt_net_sigma`). The step's source and
  the observed dH/dt stay out: the Gauss-Newton Hessian does not see the
  residual, and the source moves the stepped thickness the upwind flux
  carries by about 1e-4 on grounded ice. `--likelihood velocity` keeps the
  velocity term alone, the likelihood of the 2 km results below, whose
  spreads therefore bound the inversion's posterior from above. Every
  output records the terms and their settings, and the posterior solve
  refuses modes, and a `--resume` of saved state, of another likelihood
  (output that records none counts as the velocity term).
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
actions, nine hours at the 190 to 220 s per action the production MAP takes
on 32 ranks (2 km; the 40-mode run needed 103 of them, the 100-mode run 202
of its 404). When the restarts run out
first, the modes ARPACK has converged are kept and the json records
`converged: false`; the spreads they give are still upper bounds. `--tol`
(1e-6) is the relative accuracy asked of the eigenvalues, far finer than the
truncation.

Memory sets the other limit. At 2 km the eigensolve holds about 2.3 GB per
rank between actions and peaks near 1.4 times that during one, so 32 ranks
fit a 240 GB node; the log prints the resident set after every action, and
it must stay level. It does not with upstream tlm_adjoint 0.9.1dev, which
leaves each action's tangent-linear copies of the forward variables attached
to them, 0.6 GB per action here, which ended a 100-mode run out of memory
after 148 actions. Our tlm_adjoint (github.com/hoffmaao/tlm_adjoint, branch
`ismip7`) releases them: with it the 100-mode run's resident set rose from
62 to 66 GB over 202 actions, the growth being the eigensolver's own basis
on the root rank. `tests/test_laplace.py` fails on an environment that
still installs upstream, and both posterior scripts run the same check
(`icepack2_tools.laplace.tangent_linear_release_check`: three actions on a
small problem, counting the copies left attached) before they load the MAP,
and stop when any copy stays.

## Results at 2 km

The production MAP of the seq1 line (`map_r_friction2.h5`: friction and
fluidity controls, 1 850 366 control dofs) ran twice on one 32-rank node:
40 modes (103 Hessian actions, six hours) and 100 modes (202 actions, 12.3
hours). The leading eigenvalue is 8.4e7, the 40th 9.8e5 and the 100th 3.5e5:
every mode found so far is data-dominated by five orders of magnitude, so
the truncation still sets the spreads below, and more modes would tighten
them further. For the total grounded discharge the first 40
modes remove 77% of the prior variance and the next 60 another 14%:

| modes                  | 1   | 5   | 10  | 20  | 40  | 60  | 80  | 100 |
|------------------------|-----|-----|-----|-----|-----|-----|-----|-----|
| prior variance removed | 23% | 36% | 51% | 60% | 77% | 84% | 88% | 90% |
| sigma (Gt/yr)          | 860 | 785 | 688 | 616 | 475 | 387 | 342 | 309 |

Spreads of the grounded discharge at the MAP state, one standard deviation
in Gt/yr. The basins are the IMBIE2 basins of the ISMIP7 grid, named by
their main systems; the posterior columns are upper bounds (velocity
term alone, truncated spectrum, see above).

| quantity                            | MAP  | prior | 40 modes | 100 modes |
|-------------------------------------|------|-------|----------|-----------|
| total grounded discharge            | 1894 | 979   | 475      | 309       |
| grounding-line flux                 | 1737 | 934   | 412      | 271       |
| 0 Dronning Maud Land east           | 103  | 129   | 113      | 94        |
| 1 Enderby Land                      | 72   | 95    | 91       | 81        |
| 2 Amery                             | 70   | 138   | 88       | 74        |
| 3 Wilkes Land west (Denman)         | 127  | 201   | 201      | 69        |
| 4 Totten (Aurora)                   | 257  | 375   | 228      | 126       |
| 5 George V Land                     | 113  | 204   | 104      | 54        |
| 6 Victoria Land                     | 28   | 44    | 44       | 38        |
| 7 Ross                              | 120  | 286   | 71       | 60        |
| 8 Getz and Marie Byrd Land          | 146  | 234   | 120      | 76        |
| 9 Amundsen Sea                      | 287  | 566   | 89       | 69        |
| 10 Bellingshausen Sea               | 80   | 195   | 93       | 93        |
| 11 Peninsula west                   | 145  | 231   | 149      | 86        |
| 12 Larsen                           | 31   | 81    | 81       | 51        |
| 13 Ronne west                       | 9    | 18    | 18       | 18        |
| 14 Filchner-Ronne                   | 241  | 343   | 132      | 100       |
| 15 Dronning Maud Land west          | 64   | 105   | 70       | 67        |

Within the first 40 modes the spread falls most in the Amundsen Sea (566 to
89 Gt/yr), Ross (286 to 71), Filchner-Ronne (343 to 132), Getz (234 to 120),
the Bellingshausen Sea (195 to 93), George V Land (204 to 104), Amery (138 to
88), the western Peninsula (231 to 149) and Totten (375 to 228). Modes 41 to
100 tighten the Amundsen Sea to 69, Ross to 60, Filchner-Ronne to 100 and Getz
to 76, halve Totten (228 to 126) and George V Land (104 to 54), tighten the
western Peninsula further (149 to 86), constrain the Denman sector for the
first time (201 to 69) together with the Larsen basin (81 to 51), and add
nothing for the Bellingshausen Sea (93 at both counts). Only the small
western Ronne basin keeps its prior spread (18); Victoria Land (44 to 38)
and Enderby Land (95 to 81) move least. The prior
spread exceeds the MAP discharge in most basins, so a prior column measures
the prior's width: the linearised Gaussian prior does not know the discharge
is positive. The total's posterior spread, 309 Gt/yr or 16% of the MAP
discharge, is the number the projections start from before any of the next
steps below.

## Posterior variance by direct solve

The truncated spread falls with every mode, and the spectrum says it will go
on falling: over modes 20 to 100 the eigenvalues follow
`lambda ~ 7.8e7 i^-1.17`, which reaches 1 only near mode 5e6, beyond the
1.85 million controls. The velocity data dominate the 7.5 km prior in nearly
every direction (the production MAP's reduced chi-square is about 50 per
observed node, so a likelihood calibrated to it would still leave some 2e5
data-dominated modes), and no affordable eigensolve converges the spread.

`antarctica/scripts/run_posterior_variance.py` brackets each quantity's
variance instead. It solves `(H + Gamma^-1) x = g` by conjugate gradients
preconditioned with the prior covariance
(`icepack2_tools.laplace.posterior_variance_cg`). After `k` steps `g . x_k`
is the A-norm of the iterate, which rises monotonically to the variance: a
lower bound. Gauss-Radau quadrature on the same recurrence, with its fixed
node at 1, the bottom of the spectrum of `I + Gamma^1/2 H Gamma^1/2`, gives
an upper bound (Golub and Meurant, 2010, Sect. 6.2 and 12.2). The saved modes
seed the solve with the truncated posterior and are deflated from the
iteration, so the zero-step upper bound is the truncated variance of the
table above and every step tightens both sides. The seed costs one Hessian
action and one prior-precision action per mode; each step costs one Hessian
action at the per-action cost above, plus one precision and one covariance
action.

    ISMIP7_INVERSION=<map.h5> mpiexec -n 32 python antarctica/scripts/run_posterior_variance.py \
        --tag <tag> --quantities discharge_total,grounding_line_flux_total \
        --iterations 200 --wall-hours 22

reads the modes from `eigendec_<lc>_<tag>.h5` (`--eigendec`, `--modes`; 0
solves from the prior) and writes `posterior_<lc>_<tag>.json` after every
step: the bracket per quantity as variances and as Gt/yr, the step history,
and the file of vectors (`posterior_cg_<lc>_<tag>_<quantity>.npz`) that
`--resume` continues in a later job. `--wall-hours` stops starting steps that
long after the script began, so a 24-hour job ends on its own record. The
upper bound takes the unresolved residual to sit at the bottom of the
spectrum, so it tightens only as the solve converges; the lower bound follows
the quantity's dominant directions early. Both are bounds on the Laplace
variance under the likelihood the seed modes were computed with; the
caveats of the section above apply unchanged.

## Next steps

1. Samples of the posterior run through the forward suite, for the spread
   of the projections themselves under the production physics.
