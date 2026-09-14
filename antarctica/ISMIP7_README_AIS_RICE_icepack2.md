# ISMIP7 Projections-Antarctica README (draft)

Draft of 13 September 2026 against the June 2026 template (the Google
document linked from discussion #6). Submit one README per ice sheet, saved
as `README_AIS_RICE_icepack2`. Every item marked **[confirm]** needs a
decision or a name from the group before submission; everything else is
what the code on this branch does today.

Contributor names, affiliations and emails: **[confirm]** Andrew Hoffman
(Rice University, ah301@rice.edu) and collaborators at Indiana University and
the University of Chicago.

Date of submission: **[confirm]**

Ice Sheet Modeled/domain_id: AIS
Modeling group name/source_id: RICE **[confirm]**
Ice Sheet Model Name/ism_id: icepack2 **[confirm]** (icepack2 dual
shallow-shelf formulation on Firedrake 2026.4.1)

## Initialization methods

1. Data assimilation at a single epoch, no spin-up. The geometry
   (thickness, bed, surface) is BedMachine Antarctica sampled onto the
   model's cells; the basal friction and ice fluidity fields are inverted
   with the adjoint (L-BFGS-B, 200 iterations) as log-adjustments to a
   driving-stress Weertman anchor and a thermally derived fluidity prior.
   The misfit is the sigma-normalised velocity misfit against MEaSUReS
   (per-component errors, floored at 3 m/yr), ISSM's logarithmic velocity
   misfit, a thickness-tendency term against the observed 2003-2019 mean
   dH/dt (Smith et al. 2020, ISMIP7 observations kit) obtained from one
   implicit-Euler prognostic step with the model's own transport operator,
   and a term on the integrated net grounded mass balance. Tikhonov
   regularisation on both controls.
2. Targets: MEaSUReS v2 velocity where observed (the observation mask
   excludes the pole hole and unobserved cells), the Smith dH/dt over
   grounded ice, the integrated net mass balance. Outputs: log friction
   adjustment theta, log fluidity adjustment phi, and the velocity field
   consistent with them.
3. The velocity misfit (median 16-19 m/yr over observed nodes on the
   production mesh), the grounding-line discharge scored against the flux
   the observed velocity carries across the same facets
   (`antarctica/scripts/score_map.py`), a check that a forward re-solve
   reproduces the inversion's velocity, and a one-year balanced control
   whose volume above flotation and mass hold (`resid = 0`).
4. SMB climatology: RACMO2.4p1 (ERA5-forced) monthly SMB, 2000-2023 mean,
   used as the constant control climate and as the baseline the ISMIP7
   anomalies are added to. Ocean climatology: the ISMIP7 observational
   thermal forcing and salinity climatology (30_sep OI product). Calving
   fronts and ice margins are held at their BedMachine position during
   initialisation (there is no spin-up).
5. Yes. An apparent-mass-balance reference is computed once at the initial
   state with the model's own transport operator so that the initial
   thickness tendency is exactly zero (the ISMIP6 ctrl_proj convention of a
   balanced control); it is frozen and applied in every experiment. It is
   NOT reported in `acabf` (which is the forcing SMB) but written alongside
   it as `acabf_correction` in the same units for anyone closing the budget.

## Projections: ice-ocean and ice-shelf fracture (AIS)

6. Ocean melt: the ISMIP7 quadratic mixed-slope parameterisation of Burgard
   et al. (2022), local-quadratic variant (TF_avg = TF), with the slope
   `sin(alpha)` from the model's own draft (capped at 5e-3), constants from
   `multimelt.constants`. K is dimensionless and per IMBIE basin,
   calibrated with `antarctica/scripts/calibrate_melt.py` on the Úa 2 km
   mesh against the July 2026 re-released observation table combining
   Paolo (2023), Davison (2023) and Adusumilli (2020), integrated target
   1067.4 Gt/yr; per basin 2.4e-5 to 1.9e-4. The forward applies this
   per-basin field. K* = 4.700e-5 and the total-match K = 6.853e-5 are
   summary statistics of the fit, and the forward uses neither scalar.
   `ISMIP7_K_SCALE` multiplies the field. The previous draft quoted the
   2500 m fit against the older Paolo and Adusumilli table (865.0 Gt/yr).
   **[confirm]** that every submitted run read this calibration. Thermal
   forcing and salinity are read at the cell's draft from the ISMIP7 ocean
   forcing.
   Partially floating cells: the geometry is cell-wise (DG0); a cell is
   floating when its height above flotation is negative and then receives
   the full melt, grounded cells none. No melt on vertical ice fronts.
7. Grounding line: the flotation criterion per cell (height above
   flotation from thickness and bed); no sub-cell parameterisation.
   Basal friction is a regularised Coulomb law (`c0 = 0.5`, exact-zero on
   floating cells because the Coulomb cap is proportional to the
   effective pressure) or Budd (`N_hat = N/N_ref` with the PISM delta floor,
   gated to grounded cells by height above flotation). Effective pressure
   is the ice overburden minus the ocean pressure, so it is the migrating
   grounding line that switches friction off.
8. Calving: in the control and, by default, in the projections the calving
   front is pinned at its 2015 (BedMachine) position: ice flowing past it
   is removed and tallied as calving (`ISMIP7_FIXED_FRONT`). A level-set
   front with a von Mises calving law (Hahn, Mikula and Frolkovic 2025
   finite-volume level set; thresholds 1.0 MPa grounded, 0.15 MPa floating)
   exists but is not calibrated; **[confirm]** which the projections use.
   No sub-grid scheme beyond the sub-cell shed of the level-set front.
9. Ice-shelf collapse: the ISMIP7 collapse mask (v2.1 for CESM2-WACCM, v1
   for MRI-ESM2-0) is applied to floating cells only: every floating cell
   the year's mask flags is emptied and booked as calving
   (`ISMIP7_FRACTURE=mask`). No stress condition. **[confirm]** whether
   the projections are run with it on.
10. Tributary glaciers after a collapse: no special treatment; the front
    retreats to the new extent, the grounding line responds to the lost
    buttressing through the momentum balance, friction is unchanged.
11. The control uses the RACMO climatology and the ISMIP7 ocean climatology
    with the apparent-mass-balance correction; projections use the same
    baseline plus the ISMIP7 anomalies re-referenced to the 2000-2029 pool
    (historical 2000-2014 and ssp126 2015-2029).

## SMB questions

19. SMB is applied as a cell-mean source in the finite-volume thickness
    transport: RACMO2.4p1 climatology plus the ISMIP7 `acabf-anomaly`
    (SDBN1 8 km, v2 for CESM2-WACCM, v1 for MRI-ESM2-0) re-referenced so the
    anomaly's mean over the control window vanishes. No surface-elevation
    feedback (the `dacabfdz` gradients are not used). Forcing is NaN
    outside the downscaled mask and is filled with zero there (discussion
    #39). The last forcing year (2299 for CESM2-WACCM) is persisted for
    2300.
20. The apparent-mass-balance correction of item 5 is applied in every
    experiment, frozen at the initial state.

## GIA, bedrock and sea level

21. No bedrock adjustment.
22. BedMachine Antarctica's reference: elevations relative to the EIGEN-6C4
    geoid.
23. Not applicable.
24. Not applicable.
25. Far-field sea-level change is not included.

## Other general questions

Ice-covered area: BedMachine's mask at the initial epoch; peripheral
glaciers off the main sheet are outside the mesh; the target mask is the
initial extent, enforced at run time by the pinned front (ice past it is
removed) in the control and default projections.

PPE / ESM participation: **[confirm]**.

Summary paragraph: **[confirm, draft]** icepack2 is a finite-element
shallow-shelf model on Firedrake in its dual (velocity, membrane stress,
basal stress) formulation, with a first-order upwind finite-volume
thickness transport on the same unstructured mesh (2 km at the grounding
line, coarsening to 180 km in the interior), an adjoint initialisation to
MEaSUReS velocities and observed thickness change, regularised Coulomb
sliding, the ISMIP7 quadratic mixed-slope ocean melt with per-basin
calibration, a pinned or level-set calving front, and the ISMIP7 collapse
masks. References: Shapero et al. 2021 (icepack); Burgard et al. 2022;
Hahn, Mikula and Frolkovic 2025; Smith et al. 2020.

## Model Characteristic Table

| Characteristic | Main suite of experiments | PPE change? |
|---|---|---|
| Mesh discretisation | Delaunay triangulation (gmsh), Úa-style size field | no |
| Native grid | H: anisotropic, 2 km at the grounding line and calving front to 180 km in the interior (246,677 cells); V: vertically integrated (shallow shelf) | no |
| Native projection | EPSG:3031, same as BedMachine | no |
| Interpolation to diagnostic grid | conservative: exact cell-pixel overlap areas (supermesh) onto the 8 km grid; whole-pixel means for thickness, fluxes and fractions, covered-part means for elevations | no |
| Time integration | transport-first split: implicit Euler thickness transport, then the diagnostic solve at the new geometry; first order | no |
| Time step | 0.1 yr | no |
| Advection scheme | upwind finite volume, DG0, implicit; first order | no |
| Ice flow mechanics | shallow-shelf approximation, dual finite-element formulation (CG1 velocity, DG0 membrane and basal stress) | no |
| Ice rheology | n = 3 (composite with a linear floor for thin ice) | no |
| Basal sliding | regularised Coulomb, m = 3, c0 = 0.5 (Budd available) | **[confirm]** |
| Basal hydrology | none | no |
| Ice-shelf fracture | yes, ISMIP7 collapse mask on floating cells | no |
| Advance and retreat | grounding line free; calving front pinned at 2015 (level-set von Mises optional) | **[confirm]** |
| Grounding line | flotation criterion per cell | no |
| Calving | pinned front; ice past it removed | **[confirm]** |
| Initial SMB | RACMO2.4p1 2000-2023 climatology | no |
| Bedrock adjustment | no | no |
| Year of initial condition | 2015 | no |
| Densities, gravity | rho_i = 917, rho_o = 1024, rho_w = 1000 kg m-3; g = 9.81 m s-2 | no |
| Variables not included | none of the mandatory set; no 3D or thermal variables (no thermal model); `hfgeoubed`, `litemp*`, `zvel*`, `thdrflf`, `deltag`, `refgeoid` absent | no |
| Days per year | 365.25: the model's year is icepack's, 31557600 s, and every model-to-SI conversion in the submitted files uses it. The time axis in the files is the standard calendar | no |
| Other | apparent-mass-balance correction frozen at the initial state; forcing versions cited per file in the submission | no |

Forcing versions used (from `antarctica/scripts/audit_forcing_versions.py`,
13 September 2026): CESM2-WACCM atmosphere SDBN1-8000m v2, ocean v3,
fracture v2.1; MRI-ESM2-0 atmosphere GEMB-SDBN1-8000m v1, ocean v3,
fracture v1; ISMIP7 ocean climatology 30_sep; observations kit
AntarcticaObsISMIP7-v1.2.
