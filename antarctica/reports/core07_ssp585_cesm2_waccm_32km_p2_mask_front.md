# Core 7: ssp585_cesm2_waccm_p2_mask_front (32 km)

- date: 2026-09-22
- git: a0f92a8
- log: `../logs/ismip7_fwd_10569015.out`
- timeseries: `results/ssp585_cesm2_waccm_p2_mask_front_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly CESM2-WACCM ssp585 SDBN1-8000m v2
- Forcing provenance: atmosphere acabf CESM2-WACCM ssp585 SDBN1-8000m v2
- Forcing provenance: ocean tf CESM2-WACCM ssp585 ocean v3
- Forcing provenance: ocean so CESM2-WACCM ssp585 ocean v3
- Forcing provenance: fracture collapse_mask CESM2-WACCM ssp585 ice_shelf_collapse_mask_cesm2waccm_ssp585_ismip7_8km-v2.1.nc v2.1
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=mask_front (flagged floating cells are removed once open water reaches them, and booked as calving)
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=mask_front ended t=2301.0 with flagged=5158 removed=5154 held=4 cells; most held at once 64 cells at t=2103.1 (steps from t=2015.0)

## Run environment

```
ISMIP7_APPARENT_MB=1
ISMIP7_AUTO_RESUME=1
ISMIP7_BUFFER_M=0
ISMIP7_CHECKPOINT_EVERY_YR=5
ISMIP7_CLIM_END=2029    # default (not exported)
ISMIP7_CLIM_SCENARIO=ssp126    # default (not exported)
ISMIP7_CLIM_START=2000    # default (not exported)
ISMIP7_CONTINUATION_STEPS=8    # default (not exported)
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=full_mumps
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER_CANONICAL=full_mumps    # resolved canonical mode
ISMIP7_DT=0.1
ISMIP7_EXPERIMENT=ssp585_cesm_waccm
ISMIP7_FIXED_FRONT=1
ISMIP7_FRACTURE=mask_front
ISMIP7_FRICTION=budd
ISMIP7_GEOMETRY_SPACE=dg0
ISMIP7_KEEP_CHECKPOINTS=3
ISMIP7_KSP_MAXIT=1000    # default (not exported)
ISMIP7_KSP_RTOL=1e-6    # default (not exported)
ISMIP7_K_MELT=8.5e-05    # default (not exported)
ISMIP7_K_PER_BASIN_NPZ=/N/project/ice_rheology/ISMIP7/antarctica/results/issue11_melt_check/K_issue11_mesh2500.npz
ISMIP7_LC=32000
ISMIP7_LC_COARSE=320000
ISMIP7_MASS_RESIDUAL_TOL_GT=5e-5    # default (not exported)
ISMIP7_MELT_SLOPE=local
ISMIP7_N_FLOW=3
ISMIP7_OUTPUT=1
ISMIP7_OUTPUT_INTERVAL=10
ISMIP7_RESCUE_ENABLED=1    # default (not exported)
ISMIP7_RESCUE_MAXIT=600    # default (not exported)
ISMIP7_RUN_TAG=p2_mask_front
ISMIP7_SIN_ALPHA_ANT=0.005115    # default (not exported)
ISMIP7_SNES_ATOL=1e-50    # default (not exported)
ISMIP7_SNES_ATOL_SCALE=100    # default (not exported)
ISMIP7_SNES_DIVERGENCE_TOL=-3    # default (not exported)
ISMIP7_SNES_KSP_EW=0    # default (not exported)
ISMIP7_SNES_LINESEARCH=nleqerr    # default (not exported)
ISMIP7_SNES_LOG=stdout    # default (not exported)
ISMIP7_SNES_MAXIT=200    # default (not exported)
ISMIP7_SNES_MONITOR=0    # default (not exported)
ISMIP7_SNES_RESTART_FAILURE_ATOL_SCALE=1e-6    # default (not exported)
ISMIP7_SNES_RTOL=1e-8    # default (not exported)
ISMIP7_SNES_STOL=0    # default (not exported)
ISMIP7_SNES_TYPE=newtonls    # default (not exported)
ISMIP7_SOLVER_VIEW=0    # default (not exported)
ISMIP7_SUBCYCLES=1,4,16,64
ISMIP7_TRANSPORT_KSP_MAXIT=500    # default (not exported)
ISMIP7_TRANSPORT_KSP_RTOL=1e-10    # default (not exported)
OMP_NUM_THREADS=1
```

## Solver configuration

```json
{
  "continuation_steps": 8,
  "diagnostic_label": "full-jacobian-mumps",
  "diagnostic_mode": "full_mumps",
  "diagnostic_mode_requested": "full_mumps",
  "diagnostic_petsc_options": {
    "ksp_type": "gmres",
    "mat_mumps_cntl_3": 1e-12,
    "mat_mumps_icntl_14": 400,
    "mat_mumps_icntl_24": 1,
    "mat_type": "aij",
    "pc_factor_mat_solver_type": "mumps",
    "pc_type": "lu",
    "snes_atol": 1e-50,
    "snes_divergence_tolerance": -3.0,
    "snes_linesearch_type": "nleqerr",
    "snes_max_it": 200,
    "snes_rtol": 1e-08,
    "snes_stol": 0.0,
    "snes_type": "newtonls"
  },
  "linearization_state": "assembled",
  "mass_residual_tolerance_gt": 5e-05,
  "monitoring": {
    "destination": "stdout",
    "enabled": false,
    "view": false
  },
  "options_prefixes": {
    "diagnostic": "ismip7_diagnostic_",
    "transport": "ismip7_transport_"
  },
  "rescue_enabled": true,
  "rescue_max_it": 600,
  "snes_atol_policy": {
    "initial": 1e-50,
    "post_convergence_scale": 100.0,
    "restart_failure_residual_scale": 1e-06
  },
  "subcycles": [
    1,
    4,
    16,
    64
  ],
  "transport_petsc_options": {
    "ksp_error_if_not_converged": null,
    "ksp_max_it": 500,
    "ksp_rtol": 1e-10,
    "ksp_type": "gmres",
    "pc_type": "bjacobi",
    "sub_ksp_type": "preonly",
    "sub_pc_type": "ilu"
  }
}
```

## Budget at marker years

year | vaf_mm_sle | mass_gt | smb_gtyr | melt_gtyr | outflux_gtyr | calv_gt | clamp_gt | resid_gt | amb_gtyr | collapse_flagged_cells | collapse_removed_cells | collapse_held_cells
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
2015.1 | 56995.801924 | 23687709.37 | 2483.1482 | 1166.9617 | 1638.0124 | 0.0549 | 4.1261 | -0.0000 | 554.2480 | 0 | 0 | 0
2055.9 | 57064.811108 | 23686758.10 | 2543.6681 | 2327.6152 | 1331.6055 | 3.3840 | 11.8148 | -0.0000 | 554.2480 | 87 | 77 | 10
2096.8 | 57179.139275 | 23633168.08 | 2857.7398 | 5090.8425 | 683.8274 | 40.9327 | 157.3137 | 0.0000 | 554.2480 | 790 | 737 | 53
2137.6 | 57348.674546 | 23522737.56 | 1791.4950 | 15154.5865 | 361.9334 | 576.5154 | 1151.1726 | -0.0000 | 554.2480 | 2003 | 1998 | 5
2178.5 | 57489.883089 | 23381934.96 | -1210.4967 | 32340.0876 | 209.2378 | 92.2951 | 3166.7315 | -0.0000 | 554.2480 | 3244 | 3228 | 16
2219.3 | 57561.125727 | 23251724.95 | -4495.4009 | 59831.9882 | 202.9203 | 126.3009 | 6437.8218 | -0.0000 | 554.2480 | 4239 | 4236 | 3
2260.2 | 57551.365791 | 23163085.32 | -7745.4860 | 101511.0558 | 193.0264 | 144.9549 | 10925.2658 | 0.0000 | 554.2480 | 4799 | 4794 | 5
2301.0 | 57475.818324 | 23098512.01 | -9491.5655 | 124570.7020 | 183.9365 | 148.1089 | 13410.2590 | -0.0000 | 554.2480 | 5158 | 5154 | 4

## Observational audit

```
ISMIP6-track audit: ssp585_cesm2_waccm_p2_mask_front_32000_timeseries.csv
  2860 steps, 2015.1->2301.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                       -1373.1   [ 2000.0,  2900.0] Gt/yr     FAIL
  shelf basal melt          40559.6   [  600.0,  1800.0] Gt/yr     FAIL
  front discharge            2301.1   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)         -2068.2   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)           1.7   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway     137656.8   [yr-median < 6000, growth<1.5x for 2 yr] FAIL

  OFF TRACK (4 FAIL rows)
```
