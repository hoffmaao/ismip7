# Core 7: ssp585_cesm2_waccm_i116y2003off (32 km)

- date: 2026-09-26
- git: 9c24421
- log: `../logs/ismip7_fwd_10659477.out`
- timeseries: `results/ssp585_cesm2_waccm_i116y2003off_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly CESM2-WACCM ssp585 SDBN1-8000m v2
- Forcing provenance: atmosphere acabf CESM2-WACCM ssp585 SDBN1-8000m v2
- Forcing provenance: ocean tf CESM2-WACCM ssp585 ocean v3
- Forcing provenance: ocean so CESM2-WACCM ssp585 ocean v3
- Forcing provenance: ocean melt calibration K_issue11_mesh2500.npz sha256 1d2a66a3089bdce7c94ec825f49f66363a394f936606c4041b8c9120b3e1c32d: a legacy per-basin K named with ISMIP7_K_PER_BASIN_NPZ, not the tracked calibration
- SMB-elevation feedback: off (ISMIP7_SMB_ELEVATION_FEEDBACK=0)
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)
- Calving front owner: legacy fixed-front mask (ISMIP7_FIXED_FRONT)

## Run environment

```
ISMIP7_APPARENT_MB=1
ISMIP7_AUTO_RESUME=1
ISMIP7_BNDIDS=/N/project/ice_rheology/ISMIP7/antarctica/mesh/boundary_ids_antarctica_320000_32000_buffered0.json
ISMIP7_BUFFER_M=0
ISMIP7_CHECKPOINT_EVERY_YR=5
ISMIP7_CLIM_END=2029    # default (not exported)
ISMIP7_CLIM_SCENARIO=ssp126    # default (not exported)
ISMIP7_CLIM_START=2000    # default (not exported)
ISMIP7_CONTINUATION_STEPS=8    # default (not exported)
ISMIP7_DATA_ROOT=/N/project/ice_rheology/ISMIP7/ISMIP7/AIS
ISMIP7_DELTAT_PER_BASIN_NPZ=none, the legacy per-basin K named in ISMIP7_K_PER_BASIN_NPZ    # default (not exported)
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=full_mumps
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER_CANONICAL=full_mumps    # resolved canonical mode
ISMIP7_DT=0.1
ISMIP7_EXPERIMENT=ssp585_cesm_waccm
ISMIP7_FIXED_FRONT=1
ISMIP7_FRACTURE=none
ISMIP7_FRICTION=budd
ISMIP7_GEOMETRY_SPACE=dg0
ISMIP7_INVERSION=/N/project/ice_rheology/ISMIP7/antarctica/mesh/inversion_icepack2_budd_n3_dg0_logvelnet_32000.h5
ISMIP7_KEEP_CHECKPOINTS=3
ISMIP7_KSP_MAXIT=1000    # default (not exported)
ISMIP7_KSP_RTOL=1e-6    # default (not exported)
ISMIP7_K_PER_BASIN_NPZ=/N/project/ice_rheology/ISMIP7/antarctica/results/issue11_melt_check/K_issue11_mesh2500.npz
ISMIP7_LC=32000
ISMIP7_LC_COARSE=320000
ISMIP7_MAP_DEFAULT=/N/scratch/dlilien/ismip7_issue116/run_y2003/antarctica/mesh/inversion_icepack2_budd_n3_dg0_logvelnet_32000.h5
ISMIP7_MASS_RESIDUAL_TOL_GT=5e-5    # default (not exported)
ISMIP7_MELT_SLOPE=local
ISMIP7_MESH=/N/project/ice_rheology/ISMIP7/antarctica/mesh/antarctica_320000_32000_buffered0.msh
ISMIP7_N_FLOW=3
ISMIP7_OBS_DATA_ROOT=/N/project/ice_rheology/ISMIP7/antarctica/data
ISMIP7_OUTPUT=0
ISMIP7_OUTPUT_INTERVAL=10
ISMIP7_RESCUE_ENABLED=1    # default (not exported)
ISMIP7_RESCUE_MAXIT=600    # default (not exported)
ISMIP7_RUN_TAG=i116y2003off
ISMIP7_SIN_ALPHA_ANT=0.005115    # default (not exported)
ISMIP7_SMB_ELEVATION_FEEDBACK=0
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
2015.1 | 56845.051726 | 23647905.10 | 2482.8965 | 1350.5526 | 1713.0464 | 0.0608 | 0.0584 | 0.0000 | 630.9474 | 0 | 0 | 0
2055.9 | 56882.523447 | 23638601.14 | 2543.2337 | 2522.7451 | 1537.3766 | 0.0638 | 1.9862 | 0.0000 | 630.8955 | 0 | 0 | 0
2096.8 | 56952.467911 | 23584726.35 | 2857.2913 | 5480.8967 | 986.4541 | 0.0329 | 62.5663 | 0.0000 | 627.6086 | 0 | 0 | 0
2137.6 | 57043.983572 | 23446384.17 | 1811.9074 | 15128.8410 | 474.0721 | 0.0250 | 792.2911 | 0.0000 | 639.1219 | 0 | 0 | 0
2178.5 | 57089.761785 | 23240794.98 | -997.7199 | 29068.3186 | 239.7164 | 0.0204 | 2506.6174 | 0.0000 | 594.4609 | 0 | 0 | 0
2219.3 | 57042.315473 | 23067986.39 | -4147.3850 | 52937.4025 | 162.6456 | 0.0158 | 5318.3075 | -0.0000 | 539.0428 | 0 | 0 | 0
2260.2 | 56910.632720 | 22933006.74 | -7318.3347 | 90305.8673 | 129.2299 | 0.0110 | 9265.5103 | -0.0000 | 602.5783 | 0 | 0 | 0
2301.0 | 56721.381335 | 22787330.97 | -8960.3000 | 103285.2205 | 107.0422 | 0.0105 | 10915.8460 | 0.0000 | 612.6983 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ssp585_cesm2_waccm_i116y2003off_32000_timeseries.csv
  2860 steps, 2015.1->2301.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                       -1203.4   [ 2000.0,  2900.0] Gt/yr     FAIL
  shelf basal melt          36215.6   [  600.0,  1800.0] Gt/yr     FAIL
  front discharge             624.4   [  700.0,  2400.0] Gt/yr     WARN
  dM/dt (post-2016)         -3019.7   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)          -0.4   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1719.1   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  OFF TRACK (3 FAIL rows)
```
