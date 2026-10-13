# Staged initialization, relaxation and forward routine

A written-out sequence from BedMachine and MEaSUREs to the ISMIP7 forwards,
built from the sequential initialization several ISMIP6 groups used (shelf
rheology first, then grounded friction; Seroussi et al. 2020, Appendix C)
and a free relaxation before the historical, on this repository's own tools.
Every stage is one existing driver under environment knobs; the recipe in
`scripts/initialization_routine.sh` runs them in order and records each
stage's MAP or state. The 32 km numbers quoted are from the 13 October 2026
test matrix (`project` run records `test-32km-framework-*`).

The routine answers two defects of the joint descent: the friction gradient
is ~50x the fluidity gradient per unit control, so a joint descent fits the
friction first and leaves the shelves at their prior for tens of iterations
and then dies by line search; and a MAP fitted to velocity alone carries a
flux divergence the forward turns into drift, which the apparent mass
balance freezes rather than removes.

## Stage 0: geometry and data conditioning

| what | knob | why |
|---|---|---|
| BedMachine as hat-weighted patch means | `ISMIP7_RASTER_SAMPLE=patch` (`patch_front` on a coarse mesh) | conserves the raster volume and is smoother than the vertex sample (grounded surface jump -11 % at 2 km); the front rule is only needed at 32 km |
| MEaSUREs as patch means with the patch rms error as sigma | `ISMIP7_VEL_SAMPLE=patch` | the datum a CG1 node can represent; the point sample's representativeness error (23 m/yr median at 32 km, 1.7 at 2 km) no longer scores as model error |
| detached ice removed | `ISMIP7_KILL_BERGS=1` | a fragment carries the terminus push on every facet and nothing holds it |
| thickness floor | `ISMIP7_H_CLAMP=10` | ice-free cells that touch ice at a vertex hold 10 m instead of a cliff |
| subglacial lake base | `ISMIP7_LAKE_ICE_BASE=1` (default) | |
| speed cap | `ISMIP7_U_OBS_MAX=5000` | the fastest pixels are the noisiest; inert at 32 km, to be measured at 2 km |

Together the margin conditioning alone took the 32 km velocity chi^2 from
373 to 348 per node at 120 iterations. A patch-sampled run needs its melt
calibration fitted under patch sampling (`calibrate_deltaT.py` with
`ISMIP7_RASTER_SAMPLE=patch`, named with `ISMIP7_DELTAT_PER_BASIN_NPZ`): the
melt contract refuses a fit made under another sampling.

## Stage 1: shelf rheology, bounded

`ISMIP7_INVERT=phi ISMIP7_FLUIDITY_CONTROL=floating ISMIP7_PHI_BOUNDS=-2.3,2.3`

Only the fluidity on floating ice moves; grounded ice keeps the prior
fluidity (the depth-averaged temperature raster) and the friction stays at
its start. The bounds hold A within exp(+-2.3) of the prior, the range a
temperature between -1 and -60 C spans; without them a shelf-first stage from
a cold start ran phi away (prior term 1e6 within 13 iterations, each
evaluation three times slower). 40 iterations.

## Stage 2: grounded friction

`ISMIP7_INVERT=theta ISMIP7_FLUIDITY_CONTROL=floating`, warm start = stage 1

Only the friction control moves. The observed thinning enters as a weight
on the velocity misfit rather than as a tendency to match
(`ISMIP7_DHDT_REWEIGHT=24 ISMIP7_DHDT_WEIGHT=0`): where the Smith field
thins fastest the velocity fit counts 25x. The tendency term
(`ISMIP7_DHDT_WEIGHT=1` with `ISMIP7_DHDT_FILTER_L`) can be kept instead,
or as well; the re-weighting changes what the MAP fits, the term changes
what the forward starts from. 40 iterations.

## Stage 3: repeat

Stages 1 and 2 again at 20 iterations each, each warm-started from the
last. At 32 km the four stages (40/40/20/20) reached a velocity chi^2 of
253 per node where the joint descent had 373 after the same 120
iterations, and the dH/dt term 8.4 against 10.1.

## Stage 4: free relaxation

`control/run.py` with `ISMIP7_APPARENT_MB` unset, `ISMIP7_FIXED_FRONT=1`,
`ISMIP7_T_END` 20 to 50 years on, the RACMO climatology and the
climatology ocean with the calibrated melt; or, for a backdated year,
`relaxation/run.py` (`ISMIP7_RELAX_START`).

Thickness, grounding line and shelves evolve freely under constant observed
forcing; the front is pinned. The year to take is where the drift of VAF,
thickness and grounded area has plateaued (the budget CSV, `dM/dt` and
`dVAF/dt` by year; 20 of 50 years in the reference workflows). The state
at that year is the relaxed geometry.

Measured at 32 km on a 120-iteration MAP: 20 years free gained mass at
+779 Gt/yr (melt falling 1244 to 762 Gt/yr as the shelves thinned, outflux
825 to 722). A free relaxation does not remove the imbalance of the
unconverged MAP; it moves it.

## Stage 5: re-invert on the relaxed geometry

Stages 1 to 3 again with `ISMIP7_WARM_START` = the relaxed state's
checkpoint. The inversion takes the geometry, the velocity observations
and the mixed state from it. The relaxed MAP carries
`geometry_source_method` and the relaxation record, and a forward on the
MAP's mesh starts from the relaxed geometry (`icepack2_tools/relaxation.py`).

Measured at 32 km: the re-inversion fitted velocity better (213 against
373 per node, with twice the total iterations) and a free forward from it
drifted +899 Gt/yr against +779 from the raw MAP, with an outflux 40 %
higher. On a MAP this far from converged the relaxation did not pay; the
2 km test is the one that counts.

## Stage 6: forced relaxation year

One year under the ESM's own forcing (`experiment.py` historical from the
relaxed state, `ISMIP7_T_END` one year on, `ISMIP7_APPARENT_MB` unset),
to absorb the jump from constant climatology to time-varying forcing before
the historical starts. Its endpoint is the historical's start. In the recipe
this is the first year of the historical itself: the historical driver runs
from 2002 (`ISMIP7_ROUTINE_HIST_START`) under one run tag, so its 2015
endpoint is the file the control and the projections look for, and no
separate restart has to be threaded through the chain.

## Stage 7: historical

`experiment.py` from 2003 (or the relaxed year), the ESM's atmosphere and
ocean, the elevation feedback on the ESM's own gradient, the observed front
as a retreat-only level set (`ISMIP7_FRONT_ADVANCE=none`, `ISMIP7_CALVING=hfb`,
`ISMIP7_FRACTURE=mask_front`), `ISMIP7_APPARENT_MB` unset. Assessment:
`check_ismip6_track.py` against the IMBIE envelope and
`compare_hist_dhdt.py` against the Smith field by basin
(`tendency_by_basin.py`). No regional rescaling of the SMB to match an
observed mass change: a dynamic loss hidden inside the SMB would
contaminate the projection anomalies.

## Stage 8: control, projections, OCX

Controls branch from the historical's 2015 endpoint (`control/run.py`,
`ISMIP7_APPARENT_MB=1` for the balanced control the protocol asks for), the
six projections from the same endpoint (`experiment.py`), OCX from the
relaxed year (`projections/ocx.py`). Submission: `ismip7_output.py`, the
scalar tool, the checker (`README.md`, section on the submission chain).

## The recipe

`scripts/initialization_routine.sh OUTDIR init` runs the stages above in
order from the shell's exported objective (the same `ISMIP7_*` knobs the
inversion chain scripts export), which it records in `OUTDIR/routine.env`
so every stage, whenever it starts, minimises or integrates the same thing.
The stage products live in `OUTDIR` (`map_<stage>.h5` with a `.done`
marker, `relaxed_state.h5`, `map_forward.h5`) and the forwards' in the
checkout's `antarctica/results`. The routine is a state machine over those
files: `next` starts every stage whose inputs exist and which has not been
started, `status` lists them, `resume STAGE` restarts a stage whose chain
died. Two launch modes, chosen at `init`:

- `ISMIP7_ROUTINE_LAUNCH=mpiexec` runs each stage in the foreground on
  `ISMIP7_ROUTINE_NRANKS` ranks (the 32 km tests);
- `ISMIP7_ROUTINE_LAUNCH=slurm` submits each stage as a chain of batch
  runner links with the resources in `ISMIP7_ROUTINE_SBATCH_INV` and
  `ISMIP7_ROUTINE_SBATCH_FWD`. An inversion stage's iterations are a budget
  across its links (`ISMIP7_CHAIN_ITER_BUDGET`, kept in `<map>.iters` by
  `inversion.sbatch`), and the link that finishes a stage calls `next`
  itself (`ISMIP7_CHAIN_THEN_CMD`), so a routine on one-hour preemptible
  links runs to the last projection with no process of its own.

The knobs: `ISMIP7_ROUTINE_ITERS` (40 for the first cycle, half for the
second), `ISMIP7_ROUTINE_RELAX_YR` (20), `ISMIP7_ROUTINE_WARM_START` (the
MAP stage 1 starts from; with no inversion stage listed, the MAP the
forwards run from), `ISMIP7_ROUTINE_ESMS`, `ISMIP7_ROUTINE_SCENARIOS`,
`ISMIP7_ROUTINE_STAGES` (a subset; `reinvert`, `hist`, `control` and `proj`
expand), `ISMIP7_ROUTINE_PROJ_AMB` (the projections' apparent mass balance,
1 like the control's, so the two share the reference built at their common
2015 state). The relaxation's end state is the control driver's final
checkpoint stamped with the relaxation record of `icepack2_tools/relaxation.py`
(the source MAP, its objective settings, the window and the step), so the
re-inversion takes the path a relaxed state has: it keeps the MAP's controls,
writes `map_<source>_relax2015.h5`, and a forward on the MAP's mesh starts
from the relaxed geometry.

## What is deliberately not in the routine

- A linear sliding law (p = q = 1): it converges in a third of the
  iterations at 32 km (269 against 373) but has no Coulomb limit at the
  grounding line and is a convergence experiment only.
- Accepting unconverged momentum solves in the forward
  (`ISMIP7_ACCEPT_UNCONVERGED=1`): kept as a robustness experiment for the
  2 km Newton walls; a published mass budget needs converged solves.
- Hand-edited geometry (a thicker shelf from another product, years of
  thinning subtracted): the backdating of the historicals is the explicit,
  recorded equivalent (`ISMIP7_GEOMETRY_BACKDATE`).
- Prescribed fronts without a law: the historical uses the observed front
  as a retreat-only bound on the calving law, not in place of it.
