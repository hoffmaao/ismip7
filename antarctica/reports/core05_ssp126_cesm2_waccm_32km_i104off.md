# Core 5: ssp126_cesm2_waccm_i104off (32 km)

- date: 2026-09-24
- git: 700a846
- log: `../logs/ismip7_fwd_10609878.out`
- timeseries: `results/ssp126_cesm2_waccm_i104off_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly CESM2-WACCM ssp126 SDBN1-8000m v2
- Forcing provenance: atmosphere acabf CESM2-WACCM ssp126 SDBN1-8000m v2
- Forcing provenance: ocean tf CESM2-WACCM ssp126 ocean v3
- Forcing provenance: ocean so CESM2-WACCM ssp126 ocean v3
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)
- Calving front owner: legacy fixed-front mask (ISMIP7_FIXED_FRONT)

## Run environment

```
ISMIP7_APPARENT_MB=0
ISMIP7_AUTO_RESUME=1
ISMIP7_BNDIDS=/N/project/ice_rheology/ISMIP7/antarctica/mesh/boundary_ids_antarctica_320000_32000_buffered0.json
ISMIP7_BUFFER_M=0
ISMIP7_CHECKPOINT_EVERY_YR=5
ISMIP7_CLIM_END=2029    # default (not exported)
ISMIP7_CLIM_SCENARIO=ssp126    # default (not exported)
ISMIP7_CLIM_START=2000    # default (not exported)
ISMIP7_CONTINUATION_STEPS=8    # default (not exported)
ISMIP7_DATA_ROOT=/N/project/ice_rheology/ISMIP7/ISMIP7/AIS
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=full_mumps
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER_CANONICAL=full_mumps    # resolved canonical mode
ISMIP7_DT=0.1
ISMIP7_EXPERIMENT=ssp126_cesm_waccm
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
ISMIP7_MESH=/N/project/ice_rheology/ISMIP7/antarctica/mesh/antarctica_320000_32000_buffered0.msh
ISMIP7_N_FLOW=3
ISMIP7_OBS_DATA_ROOT=/N/project/ice_rheology/ISMIP7/antarctica/data
ISMIP7_OUTPUT=0
ISMIP7_OUTPUT_INTERVAL=10
ISMIP7_RESCUE_ENABLED=1    # default (not exported)
ISMIP7_RESCUE_MAXIT=600    # default (not exported)
ISMIP7_RESTART=/N/project/ice_rheology/ISMIP7/antarctica/results/hist_cesm2_waccm_i104off_32000_final.h5
ISMIP7_RUN_TAG=i104off
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
2015.1 | 57335.081017 | 23770751.36 | 2325.9782 | 710.7899 | 1012.3662 | 0.1489 | 0.7219 | -0.0000 | 0.0000 | 0 | 0 | 0
2055.9 | 57474.881925 | 23798117.04 | 2677.2498 | 954.8786 | 926.6150 | 0.1801 | 1.1920 | 0.0000 | 0.0000 | 0 | 0 | 0
2096.8 | 57625.237029 | 23811897.93 | 2798.7666 | 1388.1663 | 845.9629 | 0.1926 | 3.7953 | 0.0000 | 0.0000 | 0 | 0 | 0
2137.6 | 57780.141698 | 23836504.76 | 2787.4482 | 1558.8393 | 804.9987 | 0.2058 | 2.6797 | 0.0000 | 0.0000 | 0 | 0 | 0
2178.5 | 57934.490278 | 23866232.29 | 2774.2693 | 1215.6360 | 772.3154 | 0.2189 | 5.6336 | 0.0000 | 0.0000 | 0 | 0 | 0
2219.3 | 58085.926612 | 23901462.69 | 2735.2550 | 946.6009 | 727.8813 | 0.2198 | 8.6256 | -0.0000 | 0.0000 | 0 | 0 | 0
2260.2 | 58245.131380 | 23943232.17 | 3000.5902 | 976.6191 | 715.6755 | 0.2279 | 4.5685 | 0.0000 | 0.0000 | 0 | 0 | 0
2301.0 | 58398.815589 | 23984353.59 | 2866.3574 | 1443.1329 | 687.2085 | 0.2362 | 10.8010 | -0.0000 | 0.0000 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ssp126_cesm2_waccm_i104off_32000_timeseries.csv
  2860 steps, 2015.1->2301.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2703.2   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1194.8   [  600.0,  1800.0] Gt/yr     PASS
  front discharge             807.9   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           747.5   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)           3.7   [   -2.0,     2.0] mm SLE/yr WARN
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1025.1   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  OFF TRACK (1 FAIL row)
```
