# Core 9: ctrl2015_cesm2_waccm_i116y2003on (32 km)

- date: 2026-09-26
- git: c6244e6
- log: `../logs/ismip7_fwd_10664216.out`
- timeseries: `results/ctrl2015_cesm2_waccm_i116y2003on_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- SMB climatology pool: none reported in `../logs/ismip7_fwd_10664216.out` (expected for a CTRL running on the RACMO climatology, which builds no ESM pool)
- Forcing provenance: atmosphere RACMO2.4p1 SMB climatology 2000-2029
- Forcing provenance: atmosphere dacabfdz CESM2-WACCM ctrl SDBN1-8000m v2
- Forcing provenance: ocean tf CESM2-WACCM ctrl ocean v3
- Forcing provenance: ocean so CESM2-WACCM ctrl ocean v3
- Forcing provenance: ocean melt calibration K_issue11_mesh2500.npz sha256 1d2a66a3089bdce7c94ec825f49f66363a394f936606c4041b8c9120b3e1c32d: a legacy per-basin K named with ISMIP7_K_PER_BASIN_NPZ, not the tracked calibration
- SMB-elevation feedback: dacabfdz from CESM2-WACCM ctrl SDBN1-8000m v2, surface change from the chain's initial state (H_init)
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
ISMIP7_ESM=CESM2-WACCM
ISMIP7_EXPERIMENT=control
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
2015.100000 | 56845.044678 | 23647902.69 | 2457.3672 | 1343.5549 | 1712.9309 | 0.0608 | 0.0580 | 0.0000 | 630.9474 | 0 | 0 | 0
2055.900000 | 56865.941755 | 23650747.37 | 2456.7320 | 1300.7403 | 1706.4337 | 0.0606 | 0.5825 | -0.0000 | 630.9401 | 0 | 0 | 0
2096.800000 | 56885.335593 | 23652882.46 | 2456.4707 | 1339.8031 | 1719.9837 | 0.0437 | 1.0354 | -0.0000 | 631.2630 | 0 | 0 | 0
2137.600000 | 56905.102292 | 23656008.72 | 2456.1538 | 1277.9432 | 1742.7258 | 0.1275 | 1.7584 | 0.0000 | 631.3872 | 0 | 0 | 0
2178.500000 | 56923.032525 | 23659082.82 | 2455.7574 | 1269.6274 | 1768.2101 | 0.0861 | 1.9420 | 0.0000 | 631.5802 | 0 | 0 | 0
2219.300000 | 56939.509100 | 23661602.12 | 2455.4156 | 1275.4389 | 1768.1548 | 0.0842 | 2.1758 | -0.0000 | 631.2741 | 0 | 0 | 0
2260.200000 | 56955.782308 | 23664515.57 | 2454.9752 | 1267.2421 | 1760.7811 | 0.0823 | 2.1652 | 0.0000 | 630.8928 | 0 | 0 | 0
2301.000000 | 56971.259855 | 23667909.34 | 2454.6378 | 1255.5183 | 1761.3265 | 0.0832 | 2.4075 | -0.0000 | 630.6822 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ctrl2015_cesm2_waccm_i116y2003on_32000_timeseries.csv
  2860 steps, 2015.1->2301.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2455.9   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1287.7   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1744.9   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)            70.1   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)           0.4   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1796.7   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  ON TRACK (0 FAIL rows)
```
