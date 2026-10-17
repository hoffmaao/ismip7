# Core 1: hist_cesm2_waccm_p2 (32 km)

- date: 2026-09-23
- git: a0f92a8
- log: `../logs/ismip7_fwd_10568488.out`
- timeseries: `results/hist_cesm2_waccm_p2_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly CESM2-WACCM historical SDBN1-8000m v2
- Forcing provenance: atmosphere acabf CESM2-WACCM historical SDBN1-8000m v2
- Forcing provenance: ocean tf CESM2-WACCM historical ocean v3
- Forcing provenance: ocean so CESM2-WACCM historical ocean v3
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)

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
ISMIP7_EXPERIMENT=hist_cesm_waccm
ISMIP7_FIXED_FRONT=1
ISMIP7_FRACTURE=none    # default (not exported)
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
ISMIP7_OUTPUT_INTERVAL=10
ISMIP7_RESCUE_ENABLED=1    # default (not exported)
ISMIP7_RESCUE_MAXIT=600    # default (not exported)
ISMIP7_RUN_TAG=p2
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
1850.1 | 56839.890603 | 23647063.02 | 2082.3264 | 889.2296 | 1746.5061 | 0.7597 | 0.0001 | -0.0000 | 554.2480 | 0 | 0 | 0
1873.7 | 56848.091404 | 23649623.76 | 2264.4397 | 916.1847 | 1700.6013 | 0.0889 | 1.0393 | 0.0000 | 554.2480 | 0 | 0 | 0
1897.2 | 56861.721891 | 23653637.70 | 2298.2211 | 986.2229 | 1743.3723 | 0.1079 | 0.5246 | 0.0000 | 554.2480 | 0 | 0 | 0
1920.8 | 56880.896188 | 23658847.12 | 2221.3743 | 818.8417 | 1723.9732 | 0.1008 | 1.5348 | -0.0000 | 554.2480 | 0 | 0 | 0
1944.3 | 56904.155521 | 23665335.71 | 2162.7859 | 812.3888 | 1666.1242 | 0.0501 | 2.3858 | -0.0000 | 554.2480 | 0 | 0 | 0
1967.9 | 56934.597766 | 23672823.90 | 2317.5378 | 978.2695 | 1660.9100 | 0.0559 | 2.7618 | -0.0000 | 554.2480 | 0 | 0 | 0
1991.4 | 56963.166242 | 23680199.68 | 2307.4866 | 937.6307 | 1658.3488 | 0.0599 | 3.6430 | 0.0000 | 554.2480 | 0 | 0 | 0
2015.0 | 56995.647426 | 23687682.05 | 2597.3649 | 1065.7731 | 1637.9839 | 0.0559 | 4.2791 | -0.0000 | 554.2480 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: hist_cesm2_waccm_p2_32000_timeseries.csv
  1650 steps, 1850.1->2015.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2286.7   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt            917.6   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1697.1   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           247.7   [ -400.0,   200.0] Gt/yr     WARN
  dVAF/dt (post-2016)           0.9   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1754.1   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  ON TRACK (0 FAIL rows)
```
