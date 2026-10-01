# Core 10: ctrl2015_mri_esm2_0_i116y2003on (32 km)

- date: 2026-09-26
- git: c6244e6
- log: `../logs/ismip7_fwd_10664430.out`
- timeseries: `results/ctrl2015_mri_esm2_0_i116y2003on_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- SMB climatology pool: none reported in `../logs/ismip7_fwd_10664430.out` (expected for a CTRL running on the RACMO climatology, which builds no ESM pool)
- Forcing provenance: atmosphere RACMO2.4p1 SMB climatology 2000-2029
- Forcing provenance: atmosphere dacabfdz MRI-ESM2-0 ctrl GEMB-SDBN1-8000m v2
- Forcing provenance: ocean tf MRI-ESM2-0 ctrl ocean v3
- Forcing provenance: ocean so MRI-ESM2-0 ctrl ocean v3
- Forcing provenance: ocean melt calibration K_issue11_mesh2500.npz sha256 1d2a66a3089bdce7c94ec825f49f66363a394f936606c4041b8c9120b3e1c32d: a legacy per-basin K named with ISMIP7_K_PER_BASIN_NPZ, not the tracked calibration
- SMB-elevation feedback: dacabfdz from MRI-ESM2-0 ctrl GEMB-SDBN1-8000m v2, surface change from the chain's initial state (H_init)
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
ISMIP7_ESM=MRI-ESM2-0
ISMIP7_EXPERIMENT=control
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
ISMIP7_MAP_DEFAULT=/N/scratch/dlilien/ismip7_issue116/run_ctrl/antarctica/mesh/inversion_icepack2_budd_n3_dg0_logvelnet_32000.h5
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
2015.100000 | 56840.316756 | 23646392.69 | 2457.3934 | 1240.6765 | 1697.8452 | 0.0918 | 0.0015 | 0.0000 | 491.8790 | 0 | 0 | 0
2055.900000 | 56847.586409 | 23647239.57 | 2456.9741 | 1244.0078 | 1662.8407 | 0.0833 | 0.3054 | -0.0000 | 491.8523 | 0 | 0 | 0
2096.800000 | 56861.755153 | 23651358.80 | 2456.1081 | 1165.1010 | 1656.7809 | 0.0642 | 0.9306 | -0.0000 | 491.8913 | 0 | 0 | 0
2137.600000 | 56879.991033 | 23656844.72 | 2454.9995 | 1139.9199 | 1694.5330 | 0.0738 | 1.4381 | 0.0000 | 492.1748 | 0 | 0 | 0
2178.500000 | 56897.913461 | 23662129.70 | 2454.1259 | 1102.5542 | 1727.5562 | 0.0758 | 1.6751 | 0.0000 | 492.5149 | 0 | 0 | 0
2219.300000 | 56913.651306 | 23667276.79 | 2453.4086 | 1105.9668 | 1742.1560 | 0.0713 | 1.6707 | -0.0000 | 492.2818 | 0 | 0 | 0
2260.200000 | 56929.096054 | 23671708.13 | 2452.7425 | 1100.3868 | 1756.4365 | 0.0633 | 1.8748 | -0.0000 | 492.5238 | 0 | 0 | 0
2301.000000 | 56944.766337 | 23676422.26 | 2452.1804 | 1087.3660 | 1757.9581 | 0.0608 | 2.1143 | 0.0000 | 492.0919 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ctrl2015_mri_esm2_0_i116y2003on_32000_timeseries.csv
  2860 steps, 2015.1->2301.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2454.7   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1146.9   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1707.9   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           105.4   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)           0.4   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1759.2   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  ON TRACK (0 FAIL rows)
```
