# A full rehearsal of the ISMIP7 matrix at 25 km (issue #138)

IU Quartz, 27 September 2026. Every stage from mesh to submission files ran on
the code the production matrix will use, at a resolution cheap enough to run
the whole set in a day: the mesh build, the melt refit, a warm-started Budd
re-inversion, the MAP checks, all eleven cores, the track audit, the core
reports, the ISMIP7 writer, the compliance checker and the scalar tool. The
matrix ran twice, without and with the apparent mass-balance reference
(issue #104): attempt A first, and attempt B in parallel once A's MRI-ESM2-0
control slowed to a rescue at every step.

Everything the rehearsal wrote is kept on Quartz project storage under
`antarctica/results/rehearsal_25km/`, with `MANIFEST.md` and `SHA256SUMS`, as
a candidate test submission: the `AIS/RICE/icepack2/CORE/C001` to `C011` tree
of each attempt with `params.nc`, the checker and scalar output, the movies,
and the native annual files, timeseries, final checkpoints and logs that can
regenerate the tree.

## Configuration

| item | value |
|---|---|
| code | `b554238`, the head of PR 131: main `4984512` with the dacabfdz SMB-elevation feedback on |
| mesh | `antarctica_250000_25000_buffered20000`, IU's build: 4,509 vertices, 7,615 cells, md5 `3e8b44b0`; its boundary sidecar is byte-identical to the production mesh's (md5 `5347f6c8`) |
| MAP | Budd, warm-started across meshes from Rice's `..._snap20260924_0241.h5`, 138 iterations; sha256 `8d11bc2f` |
| melt | K 6.5e-5 (K50, issue 26) with per-basin offsets refitted on this mesh; sha256 `384f8c5c` |
| forward | dt 0.025, scpc_gamg on 8 ranks, fixed front, no collapse (`ISMIP7_FRACTURE=none`), ISMIP7 output on, SMB-elevation feedback on, 2003 start |
| attempt A | `ISMIP7_APPARENT_MB=0`, tag `rehnoamb` |
| attempt B | `ISMIP7_APPARENT_MB=1`, tag `rehamb` |

Every submission named its knobs explicitly (`submit_r25.sh`, appendix):
`site_env.sh` defaults to the 1 km mesh and to regularized Coulomb, and
`simulation.py` names every output after `ISMIP7_LC`.

## Stages and gates

| stage | job | verdict |
|---|---|---|
| mesh (G1) | 10669600, 10669721 | built in 3 min. Cells in the builder's 5 km grounding-line band have a median longest edge of 32.9 km (p90 50.0 km), the usual 1.3 of a 25 km target. 26 percent of the 21,918 km exterior boundary borders ice, 1,069 of 1,401 boundary cells are ice free, and 767 cells take the ocean drag at t=0, so the 20 km ring holds open water along most of the coast |
| melt refit (G2) | 10669722, 10669724 | every basin roots between -0.90 and +0.66 K for 1067.38 Gt/yr; the K50 alone melts 1436.6 Gt/yr on this mesh. The forward's own callback reproduces every basin at ratio 1.0000 (basin 14 at 0.9999) |
| re-inversion (G3) | 10669726, 10669739 | the PR 123 ramp converged on its first rung; 2,031 of 4,509 target vertices lie outside the 2 km source mesh (the buffer ring and the coast) and take the stated fill. The log-velocity weight 85380.4 was held from the warm start across meshes (a fresh derivation gives 3.68e4). Objective 3.511e5 to 1.078e5 over 50 iterations; iterations 40 to 50 still lowered it 1.18 percent, so the budget was extended as planned and converged on the relative decrease at 88 more, 1.051e5. 18 min of wall time on 8 ranks |
| discharge score (G4a) | 10669747 | 2803 Gt/yr across the grounding line against 1848 with the observed velocity, ratio 1.52: 5.61 where the observed speed is under 100 m/yr, 0.77 at 100 to 500, 0.38 at 500 to 1500 and 0.32 above 1500. Recorded: 25 km cells do not carry the narrow fast outlets, and the inversion had converged |
| census and self-consistency (G4b) | 10669748 | 0 of 2,884 floating cells carry friction under the HAF gate (the old sign test would have put it on 396); the forward re-solves the MAP's velocity to rel L2 5.9e-8 |
| probe A (G4c) | 10669749 | ten steps without the reference: every direct solve converged, resid 0.00, speed max 1.65e4 m/yr falling, dM/dt -453 Gt/yr at step 1 |
| probe B (G4c) | 10708948 | ten steps with the reference: every direct solve converged, a_ref in -438.9 to +609.4 m/yr (net +300.1 Gt/yr), the largest thickness change 1.0 m falling to 0.2 m, dM/dt -143 then 0 Gt/yr |

## The matrix

Every core of both attempts reached its end year with the budget closed to
resid 0.00 on every row, and every core's files pass isschecker 0.5.1 with 0
errors (core 11 through `isschecker_ocx.py`) and the scalar comparison with
exit 0. The track audit judges against present-day envelopes, so a
projection fails its forced rows by design; the rows that gate a rehearsal,
the budget residual and the discharge runaway, pass on every core. Wall
time is on 8 ranks. Peak speed is the largest cell speed at any step.

Attempt B, with the reference, the candidate test submission (`rehamb`):

| core | experiment | wall (min) | rescued steps | VAF change (mm SLE) | mass change (Gt) | peak speed (m/yr, year) | track audit | checker errors | scalars |
|---|---|---|---|---|---|---|---|---|---|
| C001 | historical CESM2-WACCM | 5 | 0 | +2.9 | -763 | 1.8e+04 (2015) | ON TRACK | 0 | exit 0 |
| C002 | historical MRI-ESM2-0 | 6 | 0 | -0.7 | -1,122 | 1.7e+04 (2003) | ON TRACK | 0 | exit 0 |
| C003 | ssp370 CESM2-WACCM | 34 | 5 | +66.6 | -119,874 | 2.3e+04 (2071) | OFF TRACK (shelf basal melt; dM/dt (post-2016)) | 0 | exit 0 |
| C004 | ssp370 MRI-ESM2-0 | 28 | 1 | +58.4 | -10,033 | 1.9e+04 (2032) | ON TRACK | 0 | exit 0 |
| C005 | ssp126 CESM2-WACCM, resumed under scpc_mumps | 251 | 105 | +236.6 | -144,458 | 2.2e+04 (2069) | ON TRACK | 0 | exit 0 |
| C006 | ssp126 MRI-ESM2-0 | 116 | 5 | +203.5 | +43,338 | 2e+04 (2183) | ON TRACK | 0 | exit 0 |
| C007 | ssp585 CESM2-WACCM | 208 | 63 | -236.0 | -980,617 | 4.1e+06 (2280) | OFF TRACK (SMB; shelf basal melt; dM/dt (post-2016)) | 0 | exit 0 |
| C008 | ssp585 MRI-ESM2-0 | 130 | 17 | +180.8 | -747,590 | 7.9e+04 (2265) | OFF TRACK (SMB; shelf basal melt; dM/dt (post-2016)) | 0 | exit 0 |
| C009 | control CESM2-WACCM | 81 | 0 | +127.4 | +17,083 | 2.1e+04 (2025) | ON TRACK | 0 | exit 0 |
| C010 | control MRI-ESM2-0 | 77 | 4 | +99.7 | +26,970 | 2e+04 (2184) | ON TRACK | 0 | exit 0 |
| C011 | OCX | 8 | 0 | +11.0 | +4,385 | 1.7e+04 (2011) | ON TRACK | 0 | exit 0 |

Attempt A, without the reference (`rehnoamb`):

| core | experiment | wall (min) | rescued steps | VAF change (mm SLE) | mass change (Gt) | peak speed (m/yr, year) | track audit | checker errors | scalars |
|---|---|---|---|---|---|---|---|---|---|
| C001 | historical CESM2-WACCM | 8 | 1 | +26.3 | -4,368 | 1.6e+04 (2003) | ON TRACK | 0 | exit 0 |
| C002 | historical MRI-ESM2-0 | 10 | 3 | +25.6 | -3,718 | 1.6e+04 (2003) | ON TRACK | 0 | exit 0 |
| C003 | ssp370 CESM2-WACCM | 44 | 0 | +316.9 | -96,234 | 1.7e+04 (2029) | OFF TRACK (dM/dt (post-2016)) | 0 | exit 0 |
| C004 | ssp370 MRI-ESM2-0 | 52 | 0 | +323.5 | +19,096 | 1.6e+04 (2029) | ON TRACK | 0 | exit 0 |
| C005 | ssp126 CESM2-WACCM | 133 | 0 | +1009.6 | +45,190 | 2.8e+05 (2247) | ON TRACK | 0 | exit 0 |
| C006 | ssp126 MRI-ESM2-0 | 226 | 85 | +984.1 | +190,345 | 1.8e+04 (2185) | OFF TRACK (dM/dt (post-2016)) | 0 | exit 0 |
| C007 | ssp585 CESM2-WACCM | 421 | 10 | +500.9 | -681,907 | 7.9e+07 (2299) | OFF TRACK (SMB; shelf basal melt; dM/dt (post-2016)) | 0 | exit 0 |
| C008 | ssp585 MRI-ESM2-0, resumed under scpc_mumps | 360 | 17 | +1073.7 | -410,139 | 2.7e+06 (2236) | OFF TRACK (SMB; shelf basal melt; dM/dt (post-2016)) | 0 | exit 0 |
| C009 | control CESM2-WACCM | 147 | 4 | +868.6 | +155,196 | 3.1e+04 (2222) | OFF TRACK (dM/dt (post-2016)) | 0 | exit 0 |
| C010 | control MRI-ESM2-0 | 257 | 118 | +870.4 | +169,784 | 3.1e+04 (2218) | OFF TRACK (dM/dt (post-2016)) | 0 | exit 0 |
| C011 | OCX | 11 | 3 | +67.7 | +247 | 1.6e+04 (2003) | ON TRACK | 0 | exit 0 |

## Findings

1. **The OCX ocean on Quartz was a version behind.** `audit_forcing_versions.py`
   found all four members (cold, main, vary, warm) at v1 on disk against v2 on
   the Source Cooperative mirror. The first `check_melt_bound.py --ocx` read
   v1 and flagged basin 5 at 0.74 of the climatology and one 256 km block at
   0.46. `download_mirror.py` fetched v2 (12 files, 28.85 GB, beside v1); the
   reader takes the highest version, so every OCX run on Quartz reads v2 from
   then on, and the check against v2 flags nothing (basin 5 at 0.99). Core 11
   of attempt A was rerun on v2 and the v1 run kept aside. The audit reported
   no other stale product a core reads (issue #19, issue #41).
2. **Without the reference the runs drift and the solver struggles** (attempt
   A, issue #104). Both historicals gain about 26 mm SLE of VAF in 12 years;
   the CESM2-WACCM control gains 869 mm SLE by 2300, as the 32 km runs without
   the reference did. The interior thickness change is a cell-scale
   checkerboard of plus and minus 10 to 100 m between neighbouring cells (the
   C009 movie). The MRI-ESM2-0 control needed 118 rescued steps, 39 of them
   in the 41 steps from 2084.07 to 2085.07, at about 43 s a step, before it
   recovered.
3. **With the reference the drift is small and the checkerboard is gone**
   (attempt B). The historicals move VAF by +2.9 and -0.7 mm SLE, the
   CESM2-WACCM control by +127 mm by 2300, and both controls pass the track
   audit. Rescues remain possible: the CESM2-WACCM ssp585 needed one every
   third step around 2079, near a fast Amundsen cell, before it recovered.
4. **scpc_gamg stalls where scpc_mumps does not.** At 2166.7 in attempt A's
   MRI-ESM2-0 ssp585 every GAMG solve of the condensed system ran to its 1000
   iteration cap; step 6068 failed with DIVERGED_LINEAR_SOLVE after 736 s.
   From about 2238 attempt B's CESM2-WACCM ssp126 spent its time in long
   Newton solves (up to 168 iterations, 51 s) and rescues, about 0.13 model
   years a minute. Both were resumed from their last periodic checkpoint
   (2165 and 2260) under scpc_mumps, the plan's fallback, and ran at 0.09 to
   0.2 s a step, against 0.3 to 0.9 s under scpc_gamg in a healthy run of
   this mesh. On a mesh this small scpc_mumps is the faster solver, as the
   README says of coarse meshes on 16 ranks or fewer; the two runs' records
   name the switch.
5. **The melt column is the demand, and the clamp returns most of it late in
   ssp585** (issue #136). Attempt A's CESM2-WACCM ssp585 melts -28,523 Gt/yr
   at 2170 while the positivity clamp returns +25,164 Gt/yr: most of the melt
   falls on cells that hold less ice than a step of it removes. The budget
   closes to resid 0.00, and the track audit's shelf basal melt row reads the
   demand.
6. **Single-cell speed spikes pass unflagged.** The ssp585 runs peak at
   7.9e7 m/yr (attempt A's CESM2-WACCM, 2299), 4.1e6 m/yr (attempt B's, 2280.4,
   near (-300, -443) km) and 2.7e6 m/yr (attempt A's MRI-ESM2-0, 2236), and
   attempt A's ssp126 at 2.8e5 m/yr (2247); each leaves the budget closed and
   is gone within the year. They reach the submission: the checker gives
   B's C007 13 range warnings, velocities to 0.0133 m/s against accepted
   bounds of 0.0004 to 0.0008 m/s, and `strbasemag` to 1.05e6 Pa. The
   production runner arms no tripwire (issue #138).
7. **The 25 km MAP is not the production MAP's physics.** Its t=0
   grounding-line discharge is 1.52 of the observed-velocity discharge (5.6 in
   slow cells, 0.32 in fast outlets), and at 25 km the ssp585 runs gain VAF
   while losing mass (attempt B's MRI-ESM2-0 ssp585: +181 mm SLE, -747,590 Gt
   by 2300) because the grounding line cannot retreat through cells this
   coarse. The rehearsal tests the pipeline, and its sea-level numbers say
   nothing about the 1 km runs.
8. **Tooling gaps** (issue #138):
   - `core_report.py` recorded a check that could not run (exit 2) as off
     track or outside the envelope; fixed on this branch.
   - Two runner notes said a warm start on another rank count or mesh is
     refused; it is read by point location. Corrected on this branch.
   - `calibrate_deltaT.py` writes no `.source.json`, so a refit carries no
     record of its mesh or raster sampling until one is written by hand.
   - `site_env.sh` defaults to regularized Coulomb and the 1 km mesh, and
     `simulation.py` names outputs after `ISMIP7_LC`: a submission that omits
     a knob is silently misnamed or misconfigured.
   - isschecker 0.5.1 cannot check core 11; PR 139's `isschecker_ocx.py` does,
     and the rehearsal's C011 passes it with 0 errors (issue 18, closed).
   - Quartz's `ffmpeg` module has no libx264, so the movies are HEVC.

## Cost

76 jobs and 358 core-hours on IU Quartz for everything above, both attempts
included. The mesh, the refit and the checks took minutes each on debug
nodes, the 138-iteration inversion 18 min on 8 ranks, and each 286-year core
1.3 to 7 h on 8 ranks at dt 0.025, the slowest ones those that spent their
time in rescues or in GAMG. The output chain took 3 to 10 min a core and the
scalar tool 1 to 5.

## Artifacts and movies

On IU Quartz under `antarctica/results/rehearsal_25km/`:

| path | holds |
|---|---|
| `MANIFEST.md`, `SHA256SUMS` | every job and verdict of the rehearsal, and a checksum of every file |
| `mesh/`, `calibration/`, `maps/`, `mapcheck/` | the mesh and its QA, the refit and its sidecar, the MAP with its 50-iteration predecessor, the MAP checks |
| `results/`, `logs/` | the forwards' timeseries, final and periodic checkpoints and annual files, and every log |
| `rehamb/`, `rehnoamb/` | per attempt: `submission/AIS/RICE/icepack2/` (CORE/C001 to C011 and `params.nc`), `checker/`, `scalars/`, `reports/`, `movies/` |
| `scripts/` | the wrapper `submit_r25.sh` and every job script and tool the rehearsal ran, with `CODE_VERSIONS.txt` |

Attempt B's tree is the candidate test submission; its ids are the writer's
defaults (issue #38). A regeneration of B's C001 from the kept annual files
matched all 31 files exactly (job 10723631). The movies show, per year, the
thickness change since the experiment's first year and the speed, with the
ice edge and the grounding line: attempt B's C007 and C009, and attempt A's
C003, C005 and C009.

## Reproducing

Every submission went through `submit_r25.sh MODE [noamb|amb]` from a scratch
clone at `b554238`, and post-processing from a checkout of `b554238` merged
with main `60438f4`. The wrapper spells out every knob, since the runner
defaults are the 1 km mesh and regularized Coulomb:

```
common   ISMIP7_LC=25000 ISMIP7_LC_COARSE=250000 ISMIP7_BUFFER_M=20000 ISMIP7_FRICTION=budd
         ISMIP7_N_FLOW=3.0 ISMIP7_GEOMETRY_SPACE=dg0 ISMIP7_BNDIDS=<the sidecar>
forward  ISMIP7_INVERSION=<the MAP> ISMIP7_MESH=checkpoint ISMIP7_DELTAT_PER_BASIN_NPZ=<the refit>
         ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_gamg ISMIP7_DT=0.025 ISMIP7_APPARENT_MB=<0 or 1>
         ISMIP7_FIXED_FRONT=1 ISMIP7_FRACTURE=none ISMIP7_SMB_ELEVATION_FEEDBACK=1
mesh     timing_meshes.script, TIMING_LCS=25000 TIMING_RATIOS=10 TIMING_BUFFER=20000
refit    calibrate_deltaT.script, DELTAT_K=6.5e-5
inverse  submit.sh inversion, ISMIP7_WARM_START=<0241> ISMIP7_PRIOR_FORM=bilaplacian
         ISMIP7_PRIOR_SIGMA_THETA=0.3 ISMIP7_PRIOR_SIGMA_PHI=0.3 ISMIP7_PRIOR_RHO=7500
         ISMIP7_GRAD_PRECOND=mass_consistent ISMIP7_LOG_VEL_WEIGHT=auto ISMIP7_MAXITER=50
matrix   submit.sh projection, 8 ranks, 32 GiB, 12 h: hist_cesm_waccm with ISMIP7_CHAIN_THEN="control
         ssp126_cesm_waccm ssp370_cesm_waccm ssp585_cesm_waccm", the same for MRI-ESM2-0, and ocx
```
