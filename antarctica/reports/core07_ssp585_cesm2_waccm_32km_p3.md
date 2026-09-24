# Core 7: ssp585_cesm2_waccm_p3 (32 km)

- date: 2026-09-23
- git: 91bf573
- log: `../logs/ismip7_fwd_10595137.out`
- timeseries: `results/ssp585_cesm2_waccm_p3_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly CESM2-WACCM ssp585 SDBN1-8000m v2
- Forcing provenance: atmosphere acabf CESM2-WACCM ssp585 SDBN1-8000m v2
- Forcing provenance: ocean tf CESM2-WACCM ssp585 ocean v3
- Forcing provenance: ocean so CESM2-WACCM ssp585 ocean v3
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)
- Calving front owner: legacy fixed-front mask (ISMIP7_FIXED_FRONT)

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
ISMIP7_FRACTURE=none
ISMIP7_FRICTION=budd
ISMIP7_GEOMETRY_SPACE=dg0    # default (not exported)
ISMIP7_INVERSION=/N/project/ice_rheology/ISMIP7/antarctica/mesh/inversion_icepack2_budd_n3_dg0_logvelnet_32000.h5
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
ISMIP7_RESTART=/N/project/ice_rheology/ISMIP7/antarctica/results/hist_cesm2_waccm_p2_32000_final.h5
ISMIP7_RUN_TAG=p3
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
2055.9 | 57064.720520 | 23688687.51 | 2543.6681 | 2313.3146 | 1411.1420 | 0.0756 | 8.9079 | 0.0000 | 554.2480 | 0 | 0 | 0
2096.8 | 57169.424984 | 23651819.67 | 2857.7398 | 5373.5816 | 795.0078 | 0.0862 | 111.2882 | 0.0000 | 554.2480 | 0 | 0 | 0
2137.6 | 57311.268877 | 23556727.51 | 1791.4950 | 15757.1011 | 465.7702 | 0.0635 | 1026.3277 | 0.0000 | 554.2480 | 0 | 0 | 0
2178.5 | 57424.531675 | 23401294.10 | -1210.4967 | 32424.4028 | 265.8529 | 0.0629 | 2945.7553 | 0.0000 | 554.2480 | 0 | 0 | 0
2219.3 | 57471.999508 | 23253312.60 | -4495.4009 | 60363.7474 | 187.5246 | 0.0451 | 6151.2817 | 0.0000 | 554.2480 | 0 | 0 | 0
2260.2 | 57451.475409 | 23158556.81 | -7745.4860 | 101477.7419 | 151.3111 | 0.0161 | 10673.8226 | 0.0000 | 554.2480 | 0 | 0 | 0
2301.0 | 57370.760610 | 23069618.04 | -9491.5655 | 123455.7116 | 145.9034 | 0.0094 | 13089.8872 | 0.0000 | 554.2480 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ssp585_cesm2_waccm_p3_32000_timeseries.csv
  2860 steps, 2015.1->2301.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                       -1373.1   [ 2000.0,  2900.0] Gt/yr     FAIL
  shelf basal melt          40703.0   [  600.0,  1800.0] Gt/yr     FAIL
  front discharge             588.2   [  700.0,  2400.0] Gt/yr     WARN
  dM/dt (post-2016)         -2169.6   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)           1.3   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1640.2   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  OFF TRACK (3 FAIL rows)
```
