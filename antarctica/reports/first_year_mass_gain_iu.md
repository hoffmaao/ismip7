# The first-year mass gain in IU's Budd forwards: what the initial state carries

IU's Budd forwards gain ice from their first step. The 2 km historical from the
relaxed controls gains 15.5 mm of VAF from 2003 to 2015 (CESM2-WACCM, no
apparent mass balance), and its control keeps gaining about 1.3 mm a year. The
same runs carry a cell-scale pattern in dH (the "checkerboard"). This report
collects the tests IU ran on 13 and 14 October 2026 to find which parts of the
initial state drive that first-year imbalance. Run records:
`antarctica/runlog/test-*-first-year-*.json` and the ones named below.

## Setup shared by every test

- **Inversion:** IU's `refit_vf_budd` objective on Budd friction (m 3), n 3,
  DG0 geometry under `ISMIP7_RASTER_SAMPLE=vertex_front`: bilaplacian prior,
  sigma 30 for both controls, rho 7.5 km, `ISMIP7_LOG_VEL_WEIGHT=85380.44865839917`,
  `ISMIP7_MISFIT_NORM=sigma`, `ISMIP7_SIGMA_U_FLOOR=3`, `ISMIP7_DHDT_WEIGHT=0`,
  `ISMIP7_FLUIDITY_CONTROL=floating`, `ISMIP7_EXACT_FRONT=2`,
  `ISMIP7_DRAG_GATE=vertex`, `ISMIP7_RC_HVISC_FLOOR=2.5`, `ISMIP7_LAKE_ICE_BASE=1`,
  `ISMIP7_TRANSFER_FILL=extend`, scpc_mumps, 150 iterations unless stated.
  Warm starts: IU's 25 km `_wsiu` MAP and IU's 2 km `_ef2_vf` MAP, on their own
  meshes (`antarctica_250000_25000_buffered20000`,
  `antarctica_5000_2000_buffered20000`).
- **Forward:** a 1-year CESM2-WACCM historical from 2003 (the driver's 12-year
  Smith dH/dt backdate), `ISMIP7_APPARENT_MB=0`, `ISMIP7_FIXED_FRONT=1`,
  `ISMIP7_FRACTURE=none`, `ISMIP7_DT=0.025`, scpc_mumps, `ISMIP7_MAP_CLIP=0`,
  SMB-elevation feedback on, melt from each mesh's `vertex_front` deltaT file
  at K 6.5e-5 ("calibrated") or one uniform K 6.5e-5 on the legacy per-basin
  path ("uniform"). `resid` is 0.00 on every row of every run.
- **Measures:**

  | measure | definition |
  |---|---|
  | cell-scale tendency, `nb` | RMS over ice grounded at the start of each cell's first-year dH/dt minus the area mean of the grounded cells sharing a vertex with it |
  | cell-scale tendency, `bin` | the same against the cell's 25 km bin mean (2 km only) |
  | year-1 VAF gain | `vaf_mm_sle` at 2004.0 minus the first timeseries row (2003.025) |
  | fluxes | `front_flux_check.py`: the source side's thickness times the outward velocity, summed over grounding-line, floating-front and grounded-front facets of the t=0 ice |

## 1. The checkerboard is the initial imbalance carried forward

On the 2 km production historical and its control (CESM2-WACCM, `bin`), the
first year's cell-scale tendency (2.12 m/yr) decays and the pattern it leaves
correlates with it for decades:

| year | cell-scale RMS of dH since 2003 (m) | 25 km bin-mean RMS (m) | cell-scale tendency (m/yr) | correlation with the first year's tendency |
|---|---|---|---|---|
| 2004 | 2.1 | 0.7 | 2.12 | 1.00 |
| 2015 | 8.0 | 5.5 | 0.40 | 0.76 |
| 2040 (control) | 12.8 | 12.7 | 0.20 | 0.60 |
| 2073 (control) | 16.1 | 21.5 | 0.14 | 0.52 |

At 25 km (OCX protocol run, `vertex` sampling, uniform K, `nb`), DG0 and CG1
geometry show the same decay; CG1 is about three times noisier:

| geometry | first-year tendency (m/yr) | 2026 tendency (m/yr) | cell-scale RMS by 2026 (m) | correlation in 2026 |
|---|---|---|---|---|
| DG0 | 1.73 | 0.38 | 13.9 | 0.81 |
| CG1 | 5.05 | 1.79 | 54.7 | 0.80 |

In both geometry spaces the pattern is the first year's imbalance, carried
forward and decaying. The CG1 run first lost all its ice at
step 1, from a defect fixed in 74d0b84 (PR 168).

## 2. Sampling of the velocity target

`ISMIP7_VELOCITY_SAMPLE=cell_mean` (this PR) gives each CG1 node the area-weighted
mean of the MEaSUREs cell means over the cells around it; `point` (the default)
is the bilinear sample at the node. Each pair re-inverts from the same MAP.

| | 25 km point | 25 km cell_mean | 2 km point | 2 km cell_mean |
|---|---|---|---|---|
| nodes without a velocity observation | 1,857 of 4,509 | 734 | 50,220 of 969,764 | 38,693 |
| misfit at 150 iterations | 3,339 | 4,589 | 1,271 | 1,960 |
| \|grad\| | 72 | 333 | 0.40 | 1.36 |
| cell-scale tendency, nb (m/yr) | 1.76 | 1.21 | 1.45 | 1.40 |
| cell-scale tendency, bin (m/yr) | | | 1.94 | 1.92 |
| year-1 VAF gain, calibrated (mm) | +1.17 | +0.71 | +0.53 | +0.06 |
| year-1 VAF gain, uniform K (mm) | +1.27 | +0.72 | | |

At 25 km the cell-mean target removes 31 % of the cell-scale tendency (30 %
under uniform K). At 2 km it removes 3 % (nb) or 0.7 % (bin): a 2 km cell holds
about 20 pixels of the 450 m product, so a node's bilinear sample is already
close to an area mean. The VAF gain falls at both resolutions; section 5 traces
where.

## 3. Sampling of the geometry

- **`ISMIP7_RASTER_SAMPLE=cell_mean` with the cell-mean velocity target, 25 km:**
  the inversion's startup continuation converged, then its forward evaluations
  diverged (six failed forward solves, residual up to 1.2e17) and the optimizer
  stopped on failed trials after about 17 iterations, with no done marker
  (residual 7.6e10, misfit 1.87e5). This matches the README's measurement of
  `cell_mean` geometry at 2 km.
- **BedMachine smoothed before sampling, 25 km:** thickness, bed and surface
  filtered by a Gaussian of sigma 12.5 km (a normalized convolution inside
  BedMachine's ice mask, and separately outside it for the bed; the mask
  unchanged), then sampled `vertex_front` as usual, through a separate obs root.
  The inversion samples the filtered raster afresh
  (`ISMIP7_WARM_START_GEOMETRY=0`). The filter moved thickness by 165 m RMS
  over the ice pixels.

  | | point | point, smoothed BedMachine | cell_mean | cell_mean, smoothed BedMachine |
  |---|---|---|---|---|
  | misfit | 3,339 | 3,382 | 4,589 | 4,876 |
  | tendency, calibrated (m/yr) | 1.76 | 1.67 | 1.21 | 1.14 |
  | tendency, uniform K (m/yr) | 1.72 | 1.67 | 1.20 | 1.13 |
  | year-1 VAF gain, calibrated (mm) | +1.17 | +1.18 | +0.71 | +0.35 |
  | cell-scale RMS of the sampled H (m) | 150.6 | 129.9 | | |
  | cell-scale RMS of the sampled s (m) | 104.6 | 97.7 | | |

  Sub-cell geometry aliasing accounts for 3 to 6 % of the tendency, and adds
  roughly to the velocity target's share (34 % together). The cell-scale
  thickness structure a 25 km mesh resolves is real topography, which the
  filter leaves in place. The smoothed geometry holds 312 mm more VAF and
  33,000 Gt less ice, and its calibrated melt falls from 998 to 897 Gt/yr, so
  it serves as a diagnostic only.

## 4. DG1 geometry and the cell-mean target, 25 km

From IU's DG1 runs (PR 169's code at 15c50dc plus this PR's knob; unforced
1-year forwards, since DG1 refuses forcing, so only the differences between
columns compare):

| | DG0 point | DG0 cell_mean | DG1 point | DG1 cell_mean |
|---|---|---|---|---|
| misfit | 3,339 | 4,589 | 4,475 | 4,990 |
| year-1 VAF change (mm) | -3.92 | -4.49 | -3.80 | -4.26 |
| tendency, nb (m/yr) | 1.74 | 1.20 | 1.73 | 1.20 |

The cell-mean target lowers VAF by 0.46 to 0.58 mm and the tendency by 31 % in
both spaces. DG1 leaves the tendency where it is and moves VAF by 0.12 to
0.23 mm, inside the spread the MAPs' convergence allows. The two do not
compound. Results: Quartz `/N/scratch/dlilien/ismip7_dg1/25km/dg1_25km_2x2.json`.

## 5. Where the 2 km cell-mean run loses its VAF gain

Cell by cell, cell_mean minus point over the forward year
(`vsdiff_2km.py`, `vsflux_2km.py`):

| cells grounded at the start, by distance from the 2003 grounding line | VAF difference (mm) |
|---|---|
| 0 to 5 km | -0.24 |
| 5 to 20 km | -0.21 |
| 20 to 50 km | -0.06 |
| beyond 50 km | +0.01 |
| cells floating at the start | +0.02 |

The difference (-0.48 mm, -185 Gt of grounded ice) is taken from the 2003.0
state, so it differs by 0.01 mm from the year-1 gains in section 2, which start
at the first timeseries row (2003.025). It sits within 20 km of the grounding
line in every basin; the cells only the cell-mean target observes
carry -0.06 mm of it. Over the forward year the grounding-line flux is 1,775
against 1,949 Gt/yr and calving 1,378 against 1,590 Gt/yr.

At t=0 on the MAPs (`front_flux_check.py`, and `ffacets_2km.py` per facet):

| 2 km, Gt/yr | observed, point target | cell-mean target | point MAP | cell_mean MAP | point MAP without Shirase | cell_mean MAP without Shirase |
|---|---|---|---|---|---|---|
| grounding line | 2,156 | 2,167 | 2,256 | 2,396 | 2,229 | 2,368 |
| floating front | 998 | 1,282 | 1,320 | 1,518 | 1,297 | 1,495 |
| grounded front | 81 | 132 | 261 | 1,277 | 109 | 175 |

- Both targets give the same grounding-line flux; the cell_mean MAP fits its
  target there less closely (+11 % against +5 %). Rignot et al. (2019):
  2,050 plus or minus 100 Gt/yr.
- At the floating front the point target reads 998 Gt/yr because point sampling
  leaves front nodes without an observation; the cell-mean target's 1,282 sits
  inside BedMachine's own front flux (1,138 to 1,297, `bm_front_flux.py`).
  Both MAPs carry more than their targets.
- The grounded-front totals are dominated by three facets of thin unobserved
  front ice at Shirase Glacier, (1,385, 1,750) km: 1,101 Gt/yr at 0.5 to
  1.9e6 m/yr in the cell_mean MAP and 151 Gt/yr at about 2e5 m/yr in the point
  MAP. IU's `_ef2_vf` warm start already peaks at 1.4e5 m/yr there.

At 25 km the same check reads 1,560 Gt/yr through the grounding line under the
point target, 1,620 (point MAP) and 1,720 (cell_mean MAP).

## 6. Convergence

The 150-iteration MAPs above stopped short of convergence, except the 2 km point
MAP (0.1 % misfit change per 50 iterations). Each pair continues in stages from
its own MAP, every stage a warm start on the same mesh and sampling, with the
same 1-year forward after each (`submit_cv.sh`, `cv_submit_all.sh`,
`cv_collect.py`). 25 km, six stages of 250 iterations:

| stage (iterations) | point misfit | point tendency nb | point VAF gain (mm) | cell_mean misfit | cell_mean tendency nb | cell_mean VAF gain (mm) |
|---|---|---|---|---|---|---|
| 0 (150) | 3,339 | 1.76 | +1.17 | 4,589 | 1.21 | +0.71 |
| 1 (400) | 2,887 | 1.70 | +1.24 | 3,224 | 1.16 | +0.62 |
| 2 (650) | 2,804 | 1.69 | +1.26 | 2,963 | 1.18 | +0.41 |
| 3 (900) | 2,762 | 1.71 | +1.27 | 2,796 | 1.15 | +0.47 |
| 4 (1,150) | 2,708 | 1.69 | +1.31 | 2,729 | 1.16 | +0.45 |
| 5 (1,400) | 2,653 | 1.70 | +1.31 | 2,675 | 1.16 | +0.49 |
| 6 (1,650) | 2,624 | 1.70 | +1.31 | 2,597 | 1.16 | +0.51 |

The tendencies settle after the first stage and keep the 32 % gap while the
misfits fall another 10 to 20 %. The VAF gap widens from 0.46 to 0.80 mm. The
stage-6 MAPs carry 1,627 (point) and 1,692 (cell_mean) Gt/yr through the
grounding line against 1,560 under the point target, and calve 955 against
1,321 Gt/yr in the forward year. The 2 km pair (three stages of 300 iterations
for cell_mean, one for point) is recorded in
`test-2km-first-year-velocity-sample-convergence`.

## 7. icepack.fit

`icepack.fit` (icepack master, 5 October 2026) fits a CG1 field to point data by
least squares weighted by each datum's 1/sigma^2, with a Hessian penalty of
length alpha, in one global solve; it needs Firedrake 2026.4.2 or later. The
cell-mean target is close to that fit with alpha 0 and a lumped gain matrix.
On synthetic rasters (40 km square, 2 km cells, 400 m pixels,
`fit_vs_cellmean.py`):

| | point | cell_mean | fit, alpha 500 m | fit, alpha 2 km |
|---|---|---|---|---|
| smooth field, RMS error | 0.3 | 8.9 | 6.7 | 6.0 |
| with pixel noise (sigma 30), RMS error | 15.0 | 9.6 | 10.9 | 15.4 |
| step 50 to 500 m/yr, nodes either side | 50 / 500 | 205 / 458 | 137 / 585 | 198 / 490 |
| step, largest over and undershoot | 0 / 0 | 0 / 0 | 87 / 25 | 17 / 11 |
| nodes with a value in a data gap | 0 of 147 | 0 | 147 (RMS error 36) | 147 (36) |

The fit matches the cell mean's noise reduction at a short alpha, rings at
jumps, and fills every gap, so as a velocity target it needs an observation
mask of its own.

## Findings

1. The cell-scale pattern is the first year's imbalance carried forward, and
   it appears under DG0, CG1 and DG1 alike.
2. At 25 km the velocity target's point sampling accounts for about a third of
   that imbalance, and convergence keeps that share. At 2 km the share is 3 %.
3. Sub-cell geometry aliasing accounts for 3 to 6 % at 25 km.
4. The cell-mean target lowers the year-1 VAF gain at both resolutions by
   raising the discharge near the grounding line and the front. At the front
   that brings the target closer to BedMachine's own flux; at the grounding
   line the cell_mean MAP overshoots a target that matches the point one.
5. Thin unobserved front vertices at Shirase Glacier run at 1e5 to 1e6 m/yr in
   both 2 km MAPs.

## Cost

About 220 core-hours on IU Quartz for the work above (the 2 km velocity-sampling
pair about 170 of it), plus the 2 km convergence chain.

## Artifacts and reproducing

Quartz `/N/project/ice_rheology/ISMIP7/antarctica/results/first_year_mass_gain_iu/`:
`scripts/` (the wrappers and analyses named above, and `cv_jobs.txt`), `json/`
(every measure), `figs/` (`vsdiff_2km.png`, `gridnoise_tests.png`), `csv/` (the
test forwards' timeseries), `maps_25km/` (every 25 km MAP of these tests) and
`logs/`. The 2 km MAPs, the forward states and the smoothed BedMachine stay on
`/N/scratch/dlilien/ismip7_vs/`, which purges files untouched for 30 days; each
is rebuilt by its wrapper. The runs used the Quartz worktree
`/N/scratch/dlilien/ismip7_rehearsal25iu/runvs` at c62bcf0, which is this PR's
first commit applied to PR 168's branch at 74d0b84; the checkerboard runs used
288382a (PR 163's head) and 74d0b84 for the CG1 rerun.
