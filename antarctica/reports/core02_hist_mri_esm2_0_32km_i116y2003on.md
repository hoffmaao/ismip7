# Core 2: hist_mri_esm2_0_i116y2003on (32 km)

- date: 2026-09-26
- git: c6244e6
- log: `../logs/ismip7_fwd_10664221.out`
- timeseries: `results/hist_mri_esm2_0_i116y2003on_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly MRI-ESM2-0 historical GEMB-SDBN1-8000m v2
- Forcing provenance: atmosphere acabf MRI-ESM2-0 historical GEMB-SDBN1-8000m v2
- Forcing provenance: ocean tf MRI-ESM2-0 historical ocean v3
- Forcing provenance: ocean so MRI-ESM2-0 historical ocean v3
- Forcing provenance: ocean melt calibration K_issue11_mesh2500.npz sha256 1d2a66a3089bdce7c94ec825f49f66363a394f936606c4041b8c9120b3e1c32d: a legacy per-basin K named with ISMIP7_K_PER_BASIN_NPZ, not the tracked calibration
- Forcing provenance: atmosphere dacabfdz MRI-ESM2-0 historical GEMB-SDBN1-8000m v2
- SMB-elevation feedback: dacabfdz from MRI-ESM2-0 historical GEMB-SDBN1-8000m v2, surface change from the chain's initial state (H_init)
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)
- Calving front owner: legacy fixed-front mask (ISMIP7_FIXED_FRONT)

## Run environment

```
ISMIP7_APPARENT_MB=1
ISMIP7_AUTO_RESUME=1
ISMIP7_BNDIDS=/N/project/ice_rheology/ISMIP7/antarctica/mesh/boundary_ids_antarctica_320000_32000_buffered0.json
ISMIP7_BUFFER_M=0
ISMIP7_CHAIN_THEN=control
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
ISMIP7_ESM=MRI-ESM2-0
ISMIP7_EXPERIMENT=hist_mri_esm2
ISMIP7_FIXED_FRONT=1
ISMIP7_FRACTURE=none
ISMIP7_FRICTION=budd
ISMIP7_GEOMETRY_SPACE=dg0    # default (not exported)
ISMIP7_INVERSION=/N/project/ice_rheology/ISMIP7/antarctica/mesh/inversion_icepack2_budd_n3_dg0_logvelnet_32000.h5
ISMIP7_KEEP_CHECKPOINTS=3
ISMIP7_KSP_MAXIT=1000    # default (not exported)
ISMIP7_KSP_RTOL=1e-6    # default (not exported)
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
ISMIP7_RUN_TAG=i116y2003on
ISMIP7_SIN_ALPHA_ANT=0.005115    # default (not exported)
ISMIP7_SMB_ELEVATION_FEEDBACK=1
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
2003.100000 | 56841.394397 | 23647605.64 | 2441.0932 | 1166.9473 | 1765.1621 | 0.6023 | -0.0000 | 0.0000 | 491.8790 | 0 | 0 | 0
2004.800000 | 56840.717426 | 23647346.44 | 2122.6850 | 1171.4326 | 1762.3222 | 0.0874 | -0.0000 | -0.0000 | 491.8790 | 0 | 0 | 0
2006.500000 | 56840.661542 | 23647319.54 | 2316.4427 | 1180.3677 | 1749.6137 | 0.0874 | 0.0000 | 0.0000 | 491.8790 | 0 | 0 | 0
2008.200000 | 56840.820193 | 23647238.62 | 2323.2394 | 1293.7516 | 1740.4918 | 0.0881 | 0.0000 | -0.0000 | 491.8790 | 0 | 0 | 0
2009.900000 | 56840.496483 | 23646956.49 | 2414.6576 | 1299.7026 | 1727.1864 | 0.0892 | 0.0000 | 0.0000 | 491.8790 | 0 | 0 | 0
2011.600000 | 56840.286908 | 23646736.69 | 2367.6492 | 1289.9600 | 1719.2281 | 0.0907 | 0.0000 | 0.0000 | 491.8790 | 0 | 0 | 0
2013.300000 | 56840.013740 | 23646450.40 | 2461.5856 | 1335.4927 | 1713.9602 | 0.0913 | 0.0016 | 0.0000 | 491.8790 | 0 | 0 | 0
2015.000000 | 56840.307198 | 23646391.71 | 2552.9493 | 1316.7652 | 1699.5809 | 0.0918 | 0.0015 | 0.0000 | 491.8790 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: hist_mri_esm2_0_i116y2003on_32000_timeseries.csv
  120 steps, 2003.1->2015.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2408.2   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1264.3   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1734.2   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)          -110.4   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)          -0.1   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1771.2   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  ON TRACK (0 FAIL rows)
```
