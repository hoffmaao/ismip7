# Core 9: ctrl2015_cesm2_waccm_p3 (32 km)

- date: 2026-09-23
- git: 91bf573
- log: `../logs/ismip7_fwd_10595136.out`
- timeseries: `results/ctrl2015_cesm2_waccm_p3_32000_timeseries.csv` (gitignored; this report is the tracked record)
- observational audit: ON TRACK
- SMB climatology pool: none reported in `../logs/ismip7_fwd_10595136.out` (expected for a CTRL running on the RACMO climatology, which builds no ESM pool)
- Forcing provenance: atmosphere RACMO2.4p1 SMB climatology 2000-2029
- Forcing provenance: ocean OI climatology tf+so, release 30_sep, constant in time
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
ISMIP7_EXPERIMENT=control
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
2015.1 | 56995.795986 | 23687717.42 | 2457.6109 | 1061.5278 | 1638.0772 | 0.0610 | 4.2042 | -0.0000 | 554.2480 | 0 | 0 | 0
2055.9 | 57054.869059 | 23703911.42 | 2457.6109 | 1035.0133 | 1596.2188 | 0.0683 | 4.7270 | 0.0000 | 554.2480 | 0 | 0 | 0
2096.8 | 57111.582072 | 23721108.84 | 2457.6109 | 1045.0963 | 1587.4491 | 0.0952 | 6.1282 | -0.0000 | 554.2480 | 0 | 0 | 0
2137.6 | 57165.142536 | 23738764.58 | 2457.6109 | 1049.3055 | 1590.9348 | 0.0953 | 6.4114 | -0.0000 | 554.2480 | 0 | 0 | 0
2178.5 | 57217.330217 | 23756056.04 | 2457.6109 | 1063.8730 | 1594.1027 | 0.0912 | 6.5694 | 0.0000 | 554.2480 | 0 | 0 | 0
2219.3 | 57270.422880 | 23773288.24 | 2457.6109 | 1069.8732 | 1587.5668 | 0.0985 | 6.6648 | -0.0000 | 554.2480 | 0 | 0 | 0
2260.2 | 57323.371134 | 23790410.61 | 2457.6109 | 1062.2813 | 1600.9444 | 0.0945 | 6.7378 | 0.0000 | 554.2480 | 0 | 0 | 0
2301.0 | 57376.015478 | 23807528.79 | 2457.6109 | 1077.8286 | 1603.0097 | 0.0912 | 7.0914 | -0.0000 | 554.2480 | 0 | 0 | 0

## Observational audit

```
ISMIP6-track audit: ctrl2015_cesm2_waccm_p3_32000_timeseries.csv
  2860 steps, 2015.1->2301.0, dt=0.1 yr

  quantity                      run   obs/ISMIP6 envelope    verdict
  SMB                        2457.6   [ 2000.0,  2900.0] Gt/yr     PASS
  shelf basal melt           1057.4   [  600.0,  1800.0] Gt/yr     PASS
  front discharge            1596.4   [  700.0,  2400.0] Gt/yr     PASS
  dM/dt (post-2016)           419.2   [ -400.0,   200.0] Gt/yr     WARN
  dVAF/dt (post-2016)           1.3   [   -2.0,     2.0] mm SLE/yr PASS
  budget residual               0.0   [   -0.5,     0.5] Gt/yr     PASS
  no discharge runaway       1639.1   [yr-median < 6000, growth<1.5x for 2 yr] PASS

  ON TRACK (0 FAIL rows)
```
