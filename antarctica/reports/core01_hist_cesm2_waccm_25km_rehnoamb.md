# Core 1: hist_cesm2_waccm_rehnoamb (25 km)

- date: 2026-09-27
- git: b554238
- log: `../logs/ismip7_fwd_10669763.out`
- timeseries: `results/hist_cesm2_waccm_rehnoamb_25000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- SMB climatology pool: COMPLETE 30/30 yr, 2000-2029 (historical+ssp126, window 2000-2029, acabf-anomaly)
- Forcing provenance: atmosphere acabf-anomaly CESM2-WACCM historical SDBN1-8000m v2
- Forcing provenance: atmosphere acabf CESM2-WACCM historical SDBN1-8000m v2
- Forcing provenance: ocean tf CESM2-WACCM historical ocean v3
- Forcing provenance: ocean so CESM2-WACCM historical ocean v3
- Forcing provenance: ocean melt calibration deltaT_per_basin_25000_K6.500e-05.npz sha256 384f8c5cbc20981f0f6228921470beec40e952ee90279cabb5d4f2c35487e548 (named with ISMIP7_DELTAT_PER_BASIN_NPZ): K 6.500e-05 (K50) with a thermal-forcing offset per basin, fitted on antarctica_250000_25000_buffered20000 (IU's build for the 25 km rehearsal, 4,509 vertices); this run's mesh is antarctica_250000_25000_buffered20000
- Forcing provenance: atmosphere dacabfdz CESM2-WACCM historical SDBN1-8000m v2
- SMB-elevation feedback: dacabfdz from CESM2-WACCM historical SDBN1-8000m v2, surface change from the chain's initial state (H_init)
- Ice-shelf collapse forcing: ISMIP7_FRACTURE=none (no collapse mask is read and no cell is removed)
- Calving front owner: legacy fixed-front mask (ISMIP7_FIXED_FRONT)

25 km rehearsal evidence (icepack/ismip7 issue 138), not a submission result.

## Run environment

```
ISMIP7_APPARENT_MB=0
ISMIP7_AUTO_RESUME=1
ISMIP7_BNDIDS=/N/project/ice_rheology/ISMIP7/antarctica/results/rehearsal_25km/mesh/boundary_ids_antarctica_250000_25000_buffered20000.json
ISMIP7_BUFFER_M=20000
ISMIP7_CHAIN_THEN=control ssp126_cesm_waccm ssp370_cesm_waccm ssp585_cesm_waccm
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
ISMIP7_ESM=CESM2-WACCM
ISMIP7_EXPERIMENT=hist_cesm_waccm
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
ISMIP7_RUN_TAG=rehnoamb
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
2003.025000 | 57302.723301 | 23877510.60 | 2434.8854 | 985.3904 | 503.4617 | 34.9756 | 0.0000 | -0.0000 | 0.0000 | 0 | 0 | 0
2004.725000 | 57303.610588 | 23877102.98 | 2524.7222 | 870.9929 | 602.6592 | 31.3197 | 0.0000 | -0.0000 | 0.0000 | 0 | 0 | 0
2006.450000 | 57306.786140 | 23876549.58 | 2593.6332 | 1165.5981 | 626.2802 | 31.0857 | 0.0207 | 0.0000 | 0.0000 | 0 | 0 | 0
2008.150000 | 57311.265822 | 23875592.91 | 2545.0302 | 1222.4429 | 640.5373 | 32.1530 | 0.0118 | -0.0000 | 0.0000 | 0 | 0 | 0
2009.875000 | 57315.654989 | 23874764.40 | 2508.7110 | 982.2091 | 633.0240 | 31.7268 | 0.0117 | 0.0000 | 0.0000 | 0 | 0 | 0
2011.575000 | 57320.057836 | 23874154.54 | 2299.2857 | 961.4693 | 628.5837 | 31.9624 | 0.0242 | 0.0000 | 0.0000 | 0 | 0 | 0
2013.300000 | 57324.226822 | 23873580.50 | 2213.7797 | 744.9487 | 626.7667 | 31.8899 | 0.0129 | -0.0000 | 0.0000 | 0 | 0 | 0
2015.000000 | 57328.976999 | 23873142.30 | 2659.5487 | 897.4864 | 622.3116 | 31.7016 | 0.0259 | 0.0000 | 0.0000 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: hist_cesm2_waccm_rehnoamb_25000_timeseries.csv
  480 steps, 2003.0->2015.0, dt=0.025 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2504.9   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt            988.4   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1894.2   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)          -372.9   [ -400.0,   200.0] Gt/yr     PASS
  dVAF/dt (post-2016)           2.4   [   -2.0,     2.0] mm SLE/yr WARN
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1982.4   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  ON TRACK (0 FAIL rows)
```
