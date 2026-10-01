# Core 8: ssp585_mri_esm2_0_i136bon (25 km)

> **SUPERSEDED.** Stopped on 28 September 2026 at t = 2155.8, with the years through 2154 written: from 2153 the forcing file so_AIS_MRI-ESM2-0_ssp585_ocean_v3_2151-2160.nc on /N/project read at 0.34 MB/s while the rest of the tree read normally, the ranks sat in I/O wait for about 25 minutes a model year, and the job was cancelled. Superseded by core08-25km-ssp585-mriesm20-i136con, the same run on the final code.
>
> See `reports/MATRIX_STATUS.md` for which results are currently valid.

- date: 2026-09-28
- git: 2ba30b8
- log: `../logs/ismip7_fwd_10737272.out`
- timeseries: `results/ssp585_mri_esm2_0_i136bon_25000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: OFF TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly MRI-ESM2-0 ssp585 GEMB-SDBN1-8000m v2
- Forcing provenance: atmosphere acabf MRI-ESM2-0 ssp585 GEMB-SDBN1-8000m v2
- Forcing provenance: ocean tf MRI-ESM2-0 ssp585 ocean v3
- Forcing provenance: ocean so MRI-ESM2-0 ssp585 ocean v3
- Forcing provenance: ocean melt calibration deltaT_per_basin_25000_K6.500e-05.npz sha256 384f8c5cbc20981f0f6228921470beec40e952ee90279cabb5d4f2c35487e548 (named with ISMIP7_DELTAT_PER_BASIN_NPZ): K 6.500e-05 (K50) with a thermal-forcing offset per basin, fitted on antarctica_250000_25000_buffered20000 (IU's build for the 25 km rehearsal, 4,509 vertices); this run's mesh is antarctica_250000_25000_buffered20000
- Forcing provenance: atmosphere dacabfdz MRI-ESM2-0 ssp585 GEMB-SDBN1-8000m v2
- SMB-elevation feedback: dacabfdz from MRI-ESM2-0 ssp585 GEMB-SDBN1-8000m v2, surface change from the chain's initial state (H_init)
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)
- Calving front owner: legacy fixed-front mask (ISMIP7_FIXED_FRONT)

25 km evidence for icepack/ismip7 issue 136 (the rehearsal of issue 138 rerun from its historical on the issue 136 branch), not a submission result.

## Run environment

```
ISMIP7_APPARENT_MB=1
ISMIP7_AUTO_RESUME=1
ISMIP7_BNDIDS=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_25km/mesh/boundary_ids_antarctica_250000_25000_buffered20000.json
ISMIP7_BUFFER_M=20000
ISMIP7_CHECKPOINT_EVERY_YR=5
ISMIP7_CLIM_END=2029    # default (not exported)
ISMIP7_CLIM_SCENARIO=ssp126    # default (not exported)
ISMIP7_CLIM_START=2000    # default (not exported)
ISMIP7_CONTINUATION_STEPS=8    # default (not exported)
ISMIP7_DATA_ROOT=/N/project/ice_rheology/ISMIP7/ISMIP7/AIS
ISMIP7_DELTAT_PER_BASIN_NPZ=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_25km/calibration/deltaT_per_basin_25000_K6.500e-05.npz
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_gamg
ISMIP7_DIAGNOSTIC_LINEAR_SOLVER_CANONICAL=scpc_gamg    # resolved canonical mode
ISMIP7_DT=0.025
ISMIP7_ESM=MRI-ESM2-0
ISMIP7_EXPERIMENT=ssp585_mri_esm2
ISMIP7_FIXED_FRONT=1
ISMIP7_FRACTURE=none
ISMIP7_FRICTION=budd
ISMIP7_GEOMETRY_SPACE=dg0
ISMIP7_INVERSION=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_25km/maps/inversion_icepack2_budd_n3_dg0_logvelnet_25000_int250000_bilap_ws0241.h5
ISMIP7_KEEP_CHECKPOINTS=3
ISMIP7_KSP_MAXIT=1000    # default (not exported)
ISMIP7_KSP_RTOL=1e-6    # default (not exported)
ISMIP7_LC=25000
ISMIP7_LC_COARSE=250000
ISMIP7_MASS_RESIDUAL_TOL_GT=5e-5    # default (not exported)
ISMIP7_MELT_SLOPE=ant    # default (not exported)
ISMIP7_MESH=checkpoint
ISMIP7_N_FLOW=3.0
ISMIP7_OBS_DATA_ROOT=/N/project/ice_rheology/ISMIP7/antarctica/data
ISMIP7_OUTPUT=1
ISMIP7_RESCUE_ENABLED=1    # default (not exported)
ISMIP7_RESCUE_MAXIT=600    # default (not exported)
ISMIP7_RESTART=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_25km/results/hist_mri_esm2_0_rehamb_25000_final.h5
ISMIP7_RUN_TAG=i136bon
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
ISMIP7_SUBCYCLES=1,4,16    # default (not exported)
ISMIP7_TRANSPORT_KSP_MAXIT=500    # default (not exported)
ISMIP7_TRANSPORT_KSP_RTOL=1e-10    # default (not exported)
OMP_NUM_THREADS=1
```

## Solver configuration

```json
{
  "continuation_steps": 8,
  "diagnostic_label": "exact-local-condensation-gamg",
  "diagnostic_mode": "scpc_gamg",
  "diagnostic_mode_requested": "scpc_gamg",
  "diagnostic_petsc_options": {
    "condensed_field_ksp_atol": 5e-07,
    "condensed_field_ksp_gmres_restart": 100,
    "condensed_field_ksp_max_it": 1000,
    "condensed_field_ksp_rtol": 1e-12,
    "condensed_field_ksp_type": "fgmres",
    "condensed_field_mat_type": "aij",
    "condensed_field_mg_coarse_ksp_max_it": 50,
    "condensed_field_mg_coarse_ksp_rtol": 0.01,
    "condensed_field_mg_coarse_ksp_type": "gmres",
    "condensed_field_mg_coarse_pc_type": "jacobi",
    "condensed_field_mg_levels_ksp_type": "chebyshev",
    "condensed_field_mg_levels_pc_type": "jacobi",
    "condensed_field_near_nullspace": "none",
    "condensed_field_pc_gamg_parallel_coarse_grid_solver": null,
    "condensed_field_pc_type": "gamg",
    "ksp_max_it": 1000,
    "ksp_rtol": 1e-06,
    "ksp_type": "fgmres",
    "mat_type": "matfree",
    "pc_python_type": "icepack2_tools.preconditioners.ISMIP7SCPC",
    "pc_sc_eliminate_fields": "1,2",
    "pc_type": "python",
    "pmat_type": "matfree",
    "snes_atol": 1e-50,
    "snes_divergence_tolerance": -3.0,
    "snes_linesearch_type": "nleqerr",
    "snes_max_it": 200,
    "snes_rtol": 1e-08,
    "snes_stol": 0.0,
    "snes_type": "newtonls"
  },
  "linearization_state": "frozen",
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
    16
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
2015.025000 | 57302.068335 | 23876392.96 | 2503.6267 | 1129.9841 | 494.2141 | 31.0131 | 0.0223 | -0.0000 | 219.2613 | 0 | 0 | 0
2035.150000 | 57303.669468 | 23874234.24 | 2471.4023 | 1003.5138 | 489.1384 | 30.8603 | 1.0247 | 0.0000 | 218.8884 | 0 | 0 | 0
2055.275000 | 57314.637166 | 23873508.54 | 2669.7434 | 1294.4324 | 482.4440 | 30.2931 | 1.4673 | 0.0000 | 227.2542 | 0 | 0 | 0
2075.400000 | 57333.708603 | 23873688.75 | 3021.0150 | 1672.6966 | 473.2923 | 29.2482 | 2.5843 | 0.0000 | 225.8318 | 0 | 0 | 0
2095.525000 | 57359.908547 | 23859717.32 | 3117.5882 | 3566.7329 | 410.7780 | 26.2447 | 4.2986 | -0.0000 | 221.4626 | 0 | 0 | 0
2115.650000 | 57392.928442 | 23818220.41 | 2425.1037 | 5465.5014 | 259.5399 | 19.4576 | 16.1308 | -0.0000 | 229.1136 | 0 | 0 | 0
2135.775000 | 57420.471871 | 23738093.70 | 1400.7854 | 7190.1366 | 122.5003 | 12.5519 | 39.2082 | 0.0000 | 238.1698 | 0 | 0 | 0
2155.900000 | 57442.712639 | 23628999.28 | 433.4805 | 9495.4580 | 61.4957 | 9.1319 | 92.2392 | -0.0000 | 246.4053 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ssp585_mri_esm2_0_i136bon_25000_timeseries.csv
  5636 steps, 2015.0->2155.9, dt=0.025 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2273.4   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           3540.4   [  600.0,  1800.0] Gt/yr     FAIL
  front discharge            1318.9   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)         -1767.4   [ -400.0,   200.0] Gt/yr     FAIL
  dVAF/dt (post-2016)           1.0   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1743.7   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  OFF TRACK (2 FAIL rows)
```
