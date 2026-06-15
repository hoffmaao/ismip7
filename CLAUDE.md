# ISMIP7 Antarctica — working notes

Antarctica ISMIP7 submission using icepack2 / Firedrake. Companion to David Lilien's
Greenland repo (https://github.com/dlilien/ISMIP7_Greenland_Icepack).

## Layout

- `icepack2_tools/` — reusable Python utilities (mesh, forcing, eikonal, grounding,
  plotting, regrid). `forcing.py` reads ISMIP7 atmosphere/ocean/fracture data.
  `coupled.py` (untracked, WIP) sketches an ice ↔ plume coupling but references a
  `PlumeModel` class that does not yet exist in this tree.
- `antarctica/scripts/simulation.py` — shared simulation engine.
  - `setup_model()` loads mesh + BedMachine + MAP inversion fields, builds the
    3-field icepack2 diagnostic solver (V × Σ × τ), and ramps n=1→3.
  - `run_simulation()` time-steps: forcing callback → diagnostic solve →
    DG0 upwind implicit-Euler thickness transport with source `accum − ocean_melt`.
- `antarctica/scripts/projections/ssp585_cesm_waccm.py` — Core Experiment 7
  (SSP5-8.5 / CESM2-WACCM, 2015–2300) driver.
- `antarctica/mesh/` — `.msh` files, boundary-id JSON, `inversion_icepack2_<lc>.h5`
  MAP checkpoints (lc ∈ {500, 2000, 2500, 4000, 8000, 32000}).
- `antarctica/data/` — BedMachine, MEaSUREs velocity.
- `ISMIP7/AIS/` — ISMIP7 forcing data tree (CESM2-WACCM atmosphere & ocean,
  meltMIP targets, parameterisations).

## Environment knobs

- `ISMIP7_LC` (default 2500), `ISMIP7_LC_COARSE` (default 64000) — mesh resolution.
- `ISMIP7_MESH`, `ISMIP7_BNDIDS` — override mesh / boundary-id paths.
- `ISMIP7_T_END`, `ISMIP7_DT`, `ISMIP7_OUTPUT_INTERVAL` — time-stepping.
- `ISMIP7_RESTART` — restart checkpoint (defaults to `hist_*_final.h5` if present).
- `ISMIP7_K_MELT` — dimensionless K in Burgard quadratic-mixed-slope melt
  (default 1.15e-4, the K50 median from `parameter_selection_quadratic_example.ipynb`).
- `ISMIP7_H_CLAMP` — thickness floor (m).
- `ISMIP7_NO_CALVING_TERMINUS` — set to disable calving-terminus boundary term.
- `ISMIP7_DATA_ROOT` — override ISMIP7 data tree (defaults to repo `ISMIP7/AIS/`).

## Ocean melt parameterization — status (May 2026)

**Authoritative ISMIP7 form** (`multimelt.melt_functions.quadratic_mixed_slope`,
Burgard et al. 2022; reproduced from
https://github.com/ClimateClara/multimelt):

```
m = γ · sin(α) · (ρ_sw/ρ_i) · (c_po/L_i)² · β_S · S0 · g/(2|f|) · TF · |TF_avg|
```

with `γ` dimensionless (the tuned `K`), `melt_factor = (ρ_sw·c_pw)/(ρ_i·L_i)` [1/K],
and `U_factor = (c_po/L_i)·β_S·g/(2|f|)·S0` [m/s/K]. Calibration is per-percentile
(5/50/95) using the `parameter_selection_quadratic_example.ipynb` notebook in
`ismip/ismip7-antarctic-ocean-forcing`.

**Current implementation** (`icepack2_tools/forcing.py:quadratic_mixed_slope`):
matches `multimelt.melt_functions.quadratic_mixed_slope` verbatim with
`TF_avg = TF` (local-quadratic variant). Constants taken from
`multimelt.constants` (`f_coriolis = 1.4e-4` s⁻¹, `beta_S = 7.86e-4` PSU⁻¹,
`c_po = 3974` J/(kg·K), `L_i = 3.34e5` J/kg, `rho_sw = 1028`, `rho_i = 917`).
`S0` is read at draft depth from ISMIP7 `so` data via `ISMIP7Ocean.get_salinity`.
`sin(α)` is computed each step from `grad(s − h)` projected to CG1
(`forcing.py:compute_sin_alpha`).

**Sanity check (Paolo/Adusumilli):** integrated observed shelf melt ≈ 865 Gt/yr
over ~1.5 M km² (area-mean 0.63 m_ice/yr), from
`ISMIP7/AIS/parameterisations/ocean/meltobs/Melt_Paolo_Err_Adusumilli_imbie2_v3.csv`.
For TF = 2 K and γ_T = 1.5e-4 m/s the lumped form predicts ~3 m/yr — in the
right ballpark for cold shelves but undertuned for warm-cavity systems
(PIG/Thwaites ≈ 15–30 m/yr).

**Calibration script** (`antarctica/scripts/calibrate_melt.py`):

- Mesh from `inversion_icepack2_<LC>.h5` (default `LC=2000`, 281k vertices, GL-refined).
- Geometry **reinterpolated from BedMachine v4.1** (`bed`, `thickness`,
  `surface`, `mask`) — uses BedMachine's authoritative `mask == 3` for the
  floating field instead of relying on the h_clamp'd thickness in the
  inversion h5.
- Forcing: OI climatology `tf` and `so` at draft, IMBIE2 basins, all from
  the ISMIP7 8 km grid via `scipy.interpolate.RegularGridInterpolator`
  (NN). xarray.interp blows up memory with this many mesh nodes.
- Slope `sin(α) = |∇(s-h)| / sqrt(1 + |∇(s-h)|²)` on CG1, capped at
  `ISMIP7_SIN_ALPHA_CAP` (default 5×10⁻³) to suppress unstructured-mesh
  noise relative to multimelt's 8 km finite-difference slope.
- Per-basin integration via lumped mass + numpy sum (16 firedrake forms
  in a loop blows up memory).
- Closed form `K* = Σ M_obs·M_1 / σ² ÷ Σ M_1² / σ²` minimizes Term-1
  weighted residuals. Compares to Burgard K5=8.5e-5, K50=1.15e-4, K95=1.70e-4.

**Calibration result (May 2026, 2500m mesh):** scalar K\* (Term-1 weighted)
= 4.26×10⁻⁵; K (total-match) = 5.34×10⁻⁵. Integrated melt at K\* = 689 Gt/yr
vs observed 865 Gt/yr.

Per-basin K\_b ranges 1.66e-5 (basin 0, AP fringe shelves) to 1.47e-4
(basin 9, Amundsen). 3 of 16 basins land in Burgard K50–K95 envelope (3,
9, 15 — warm cavities); 12 are below K5 (cold cavities, where the local
quadratic over-predicts).

Per-basin output saved to
`antarctica/results/calibrated_K_per_basin_2500.npz`. Use via
`make_forcing_callback(..., K_per_basin_npz=...)` to get the calibrated
per-basin K stamped onto the mesh.

Mesh convergence: K\* differs by 3% between the 2km and 2.5km meshes — robust.

**Status (mid-session):**
- `antarctica/scripts/control/run.py` now applies climatology TF/so +
  per-basin K via a local callback in the script. Loads OI climatology
  files once into `RegularGridInterpolator` (z,y,x) and reuses every step.
- A 10-yr verification (`ISMIP7_DT=1.0`, T_END=2025) is running in the
  background on the 2500m mesh (`pid 107307`, nice +10, log
  `/tmp/ctrl_verify.log`). First prognostic step is taking >2.5 hr under
  heavy CPU contention from concurrent inversion/lcurve runs (~9 procs
  at ~300% CPU each). Reniced to let it crawl overnight.
- Sequential dt sweep (dt=0.1 over 2yr, then dt=0.01 over 0.5yr) is
  queued behind the dt=1.0 run.

**dt sweep results (May 2026, 2500m mesh):**

| dt | steps | wall | t_end | VAF | mass |
|---|---|---|---|---|---|
| 1.0 | 10 | ~7 hr | 2025.0 | 56638.10 | 24,098,050 |
| 0.1 | 20 | ~5 hr | 2017.0 | 56638.93 | 24,075,803 |
| 0.01 | 50 | ~6.5 hr | 2015.5 | 56665.84 | 24,080,695 |

- VAF agreement dt=1.0 vs dt=0.1 at t=2017: 5 mm (0.01%) ✓
- Mass at dt=0.1 is ~16,000 Gt lower than dt=1.0 — finer dt → smaller
  per-step melt → fewer `h_clamp=10m` resurrections, so dt=0.1's mass
  is more faithful to applied shelf melt.
- dt=0.01 shows pronounced early-transient VAF growth (+63 mm in 0.5
  yr) as the diagnostic flow relaxes after the n=1→3 continuation; dt=1.0
  and dt=0.1 average over this invisibly.

**Recommendation:** dt=0.1 for projections. Captures shelf melt
faithfully without oversampling the post-continuation transient.
dt=0.25 acceptable if 10 steps/yr is too expensive for a 285-yr run.

**Open issue (mitigated May 2026):** the `h_clamp=10 m` floor used to
add mass non-conservatively. Now `setup_model()` defaults to `h_clamp=0`
and uses a **composite rheology** in the action functional. See
`COMPOSITE_RHEOLOGY.md` for the full formulation.

The 10-yr `ISMIP7_DT=0.1` CTRL with composite + `h_clamp_init=10` +
`h_clamp=0` (advection) runs to completion but shows a runaway mass
gain (+120 000 Gt over 9 yr): the buffer cells, given an artificial
10 m initial layer by `h_clamp_init`, feed ice back into the melted
shelves through the diagnostic velocity. Composite handles h→0
correctly during the run, but the *initial* state is still polluted by
the clamp.

**Fix in progress:** re-invert `θ` and `φ` with composite rheology and
`h_clamp = 0`, so the MAP fields are consistent with the true
BedMachine geometry (h=0 over the buffered ocean region). The inverted
fields can then drive forward runs that never touch a thickness clamp.
`antarctica/scripts/inversion_icepack2.py` now uses composite α=1e-2
(stronger than forward's α=1e-4 because the inversion needs many
forward solves and SNES robustness is worth the extra bias).

**Next steps (in order):**
1. Switch projection scripts (`ssp585_cesm_waccm.py` et al) to use
   `make_forcing_callback(ocean=..., K_per_basin_npz=...)` and dt=0.1.
2. Investigate the `h_clamp` mass artifact — instrument total mass
   change vs (SMB - melt - GL flux) per step and report the residual.
3. Optional refinements: non-local `TF · |TF_avg|` using IMBIE2 basin
   masks; per-cell `f_coriolis` from `(x, y) → lat` on EPSG:3031;
   smoother slope (Gaussian filter on `s − h` before grad).

## Solver status

The 3-field diagnostic + DG0 upwind transport pipeline is what runs today.
The monolithic 4-field icepack2 prognostic (V × Σ × τ × Q_dg, Bernstein-bounded
Backward-Euler) is documented in memory but not yet wired in here — the MAP
transfer from the 3-field inversion to the 4-field system is the unresolved
piece.

## Useful references

- Burgard et al. 2022, *The Cryosphere*, "An assessment of basal melt
  parameterisations for Antarctic ice shelves."
- Jourdain et al. 2020 (ISMIP6 protocol) — earlier `γ_0 · melt_factor² · TF²`
  form, calibrated `γ_0 ≈ 9620 m/yr` for MEAN_ANT local quadratic.
- ISMIP7 ocean forcing pipeline:
  https://github.com/ismip/ismip7-antarctic-ocean-forcing
