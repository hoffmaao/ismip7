# ISMIP7 Antarctica

Antarctic ISMIP7 submission on icepack2 and Firedrake. Companion to the
Greenland repository (https://github.com/dlilien/ISMIP7_Greenland_Icepack).

The pipeline is a dependency chain: data, mesh, inversion, melt calibration,
control and projections.

---

## 0. What you need, by what you want to run

Install and download the rows your column ticks. Sizes are measured.

### 0.1 Software

| | Where from | Invert | Forward run | Adapt the mesh | Submit |
|---|---|:--:|:--:|:--:|:--:|
| **Firedrake 2026.4** (brings PETSc, MUMPS, mpi4py) | firedrakeproject.org | x | x | x | x |
| **icepack2** | github.com/icepack/icepack2 | x | x | x | x |
| **icepack** (raster interpolation onto meshes) | github.com/icepack/icepack | x | x | x | x |
| **icepack_tools** (`adapt_mesh`, `levelset`, `friction`, `grounding`) | github.com/hoffmaao/icepack_tools, private | | level-set front | x | |
| **tlm_adjoint** | github.com/jrmaddison/tlm_adjoint | x | | | |
| `xarray netCDF4 scipy rasterio pyproj shapely gmsh matplotlib` (`geopandas` only to build a mesh, section 3) | pip, into the Firedrake venv | x | x | x | x |
| `earthaccess` (NSIDC), `globus-sdk` (Globus route only) | pip | x | x | x | |
| **isschecker** (`ismip7-compliance-checker`) | github.com/ismip/ISM_SimulationChecker | | | | x |

`icepack_tools` is a separate private repository. Install it editable into the
same venv: `pip install -e /path/to/icepack_tools`.
`icepack2_tools/adapt_mesh.py` and `icepack2_tools/levelset.py` wrap it; the
rest of the repository runs without it.

The compliance checker needs Python 3.11 or newer. Check `python -V`; a
Firedrake venv often carries an older one, in which case give the checker its
own venv and call it by absolute path.

### 0.2 Data

| | Size | How | Invert | Forward run | Adapt | Submit |
|---|---|---|:--:|:--:|:--:|:--:|
| BedMachine Antarctica v4.1, MEaSUREs velocity v2 | 8 GB | `scripts/download_data.py` (Earthdata login) | x | x | x | |
| RACMO2.4p1 SMB climatology | 2 GB | same script | | x | | |
| ISMIP7 observations MIPkit v1.2 (Smith dH/dt) | 9 GB | `scripts/download_mirror.py --product ismip7-ais-observations data/mipkit/`, landing at `ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc` (`ISMIP7_OBS_KIT` overrides). `scripts/download_forcing.py --calibration` stages the same v1.2 file in the same place over Globus | `ISMIP7_DHDT_WEIGHT` | | `--from-obs` | |
| ISMIP7 forcing per ESM and scenario: SMB anomaly 7.5 GB, ocean `tf` 11 GB, `so` 6.9 GB (ssp585; historical 4.3 GB) | 25 GB each | `scripts/download_mirror.py` | | x | | |
| ISMIP7 fracture (collapse mask, lake properties, excess melt) | 3 GB per scenario | same, `data/<ESM>/<scenario>/fracture/` | | `ISMIP7_FRACTURE=mask` | | |
| Ocean OI climatology and IMBIE basin numbers | 3 GB | `scripts/download_forcing.py --ocean --calibration` | | x | | |
| Whole AIS tree (all ESMs, scenarios, `ctrl`, OCX, calibration) | 313 GB | same | | | | |
| Meshes and MAP checkpoints | 15 MB, 80 MB | sections 3 and 4, or from a colleague | | x | x | |
| Per-basin melt calibration `results/calibrated_K_per_basin_<lc>.npz` | 2 to 5 MB | `scripts/calibrate_melt.py` (section 5), or from a colleague | | x | | |

Source Cooperative carries the data-freeze copy and needs no account.

```bash
# one scenario for one ESM (about 25 GB)
python antarctica/scripts/download_mirror.py \
    data/CESM2-WACCM/ssp585/SDBN1-8000m/acabf-anomaly/ \
    data/CESM2-WACCM/ssp585/ocean/tf/ data/CESM2-WACCM/ssp585/ocean/so/
# the observations MIPkit (about 9 GB)
python antarctica/scripts/download_mirror.py --product ismip7-ais-observations data/mipkit/
# whether a local tree is current
python antarctica/scripts/audit_forcing_versions.py --scenario ssp585
```

Globus remains the archive of record and `download_forcing.py` drives it
(section 2a). Use it for anything the mirror has not synced.

### 0.3 Short paths

**Forward run from someone else's MAP.** Firedrake, icepack2, icepack,
BedMachine, MEaSUReS, RACMO, one scenario of forcing, their `.msh` and
`inversion_*.h5`, the OI climatology, the IMBIE basin numbers, and
`results/calibrated_K_per_basin_<lc>.npz`. That npz comes from
`calibrate_melt.py` or a colleague; the 2500 m one is mesh independent and
serves as the fallback for any `lc`. `control/run.py` and `projections/ocx.py`
abort without it; the ssp drivers warn and substitute a scalar `K`, which
completes with the wrong melt.

**Inversion.** Add `tlm_adjoint`, and the MIPkit for the dH/dt term.

**Mesh adaptation.** Add `icepack_tools`. The observation-driven size field
(`adapt_mesh.py --from-obs`) needs the MIPkit and MEaSUReS.

**Submission.** Add the checker in its own Python 3.11+ venv, run with
`ISMIP7_OUTPUT=1`, then `scripts/write_ismip7_output.py`. The two halves can
sit on different machines. The writer reads the annual checkpoints, which for a
286-year experiment come to roughly 14 GB, so it belongs where the run is; it
needs the forward's own environment plus netCDF4, since it reaches icepack2
through `icepack2_tools.dual_friction` and dates its time axis with cftime. The checker reads the written NetCDF,
a few GB after zlib, so it can stay wherever the newer Python is. Rice NOTS
carries the writer's dependencies; the checker runs here on the workstation.

### 0.4 Accounts

| For | Account | Where |
|-----|---------|-------|
| BedMachine, MEaSUREs (NSIDC) | NASA Earthdata (free) | https://urs.earthdata.nasa.gov/users/new |
| ISMIP7 forcing over Globus (the mirror needs none) | Globus + the ISMIP7 collection | https://app.globus.org |
| `icepack_tools` | access to the private repository | ask Andrew |
| Submitting results | an upload folder from the ISMIP7 team | email ismip6 at gmail.com with your Globus id, `AIS`, group name and `ism_id` |

Globus transfers to this machine also need Globus Connect Personal running
locally and its endpoint UUID.

### 0.5 Running on a cluster

`antarctica/scripts/batch_runners/` is site neutral. One file per cluster in
`sites/` holds the venv path, module loads, partitions, account and paths, and
`submit.sh` composes the scheduler command from it. Rice NOTS and IU Quartz
ship filled in, UChicago Midway is a stub, `sites/template.sh` is the blank.

```bash
antarctica/scripts/batch_runners/submit.sh inversion ISMIP7_LC=2000 \
    ISMIP7_LC_COARSE=5000 \
    ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_5000_2000_buffered0.msh
ISMIP7_SITE=iu_quartz antarctica/scripts/batch_runners/submit.sh projection \
    ISMIP7_EXPERIMENT=ssp585_cesm_waccm ISMIP7_OUTPUT=1 --dry-run
```

Details in `antarctica/scripts/batch_runners/readme.md`.

### 0.6 Repo layout

```
antarctica/
  data/            # BedMachine, velocity, RACMO                        [gitignored]
  mesh/            # *.msh, boundary_ids_antarctica_*.json, inversion_*.h5
  results/         # checkpoints, timeseries, logs                      [gitignored]
  reports/         # tracked per-core run records + MATRIX_STATUS.md
  scripts/         # entry points (sections 3 to 6)
    batch_runners/ # scheduler job scripts + sites/<cluster>.sh
ISMIP7/AIS/        # forcing tree the runtime reads                     [gitignored]
icepack2_tools/    # this repository's library
```

`FORWARD_RUN_READINESS.md` in this directory carries the data freeze, the
forcing-version audit, the OCX and control definitions, and the submission
checklist; read it before planning the core matrix.

The tracked files under `antarctica/mesh/` are the per-mesh
`boundary_ids_antarctica_*.json` sidecars (section 3). Legacy sidecars with no
mesh stem stay untracked: `boundary_ids.json` is the shared fallback that every
mesh build overwrites, so committing one would put a single build's map on the
fallback path for every clone.

---

## 1. Observational data (`download_data.py`)

```bash
cd antarctica
python scripts/download_data.py
```

| Dataset | Product | Auth | Lands in |
|---------|---------|------|----------|
| BedMachine Antarctica v4 | NSIDC-0756 | Earthdata | `data/bedmachine/` |
| MEaSUREs Ice Velocity v2 | NSIDC-0484 | Earthdata | `data/velocity/` |
| RACMO2.4p1 SMB | Zenodo `10.5281/zenodo.14217231` | none | `data/racmo/` |

Existing files are skipped, so re-running is cheap.

---

## 2. ISMIP7 forcing

Three things:

1. **The runtime tree `ISMIP7/AIS/`**, read by `icepack2_tools/forcing.py` and
   laid out per the protocol (section 2b). `ISMIP7_DATA_ROOT` overrides the
   root.
2. **`download_mirror.py`**, the Source Cooperative route: anonymous HTTPS,
   resumable, no Globus endpoint, and the only option on a machine without one.
   It takes product-relative prefixes and lands files under `--root` in the
   versioned layout `forcing.py` expects. Its module docstring documents the
   prefixes, the resume rule and the size check.
3. **`download_forcing.py`**, the Globus route, mirroring the collection
   straight into the runtime layout: climatology, bias and calibration files
   (`--ocean`, `--calibration`), and the per-(ESM, scenario) sets
   (`--scenarios`). `python scripts/preflight.py` then reports which core
   experiments the local tree can run.

### 2a. Using `download_forcing.py`

```bash
cd antarctica
python scripts/download_forcing.py --login        # caches ~/.ismip7_globus_tokens.json
python scripts/download_forcing.py --list         # legacy climatology subtree; /ISMIP7/AIS needs the Globus web app
export GLOBUS_LOCAL_ENDPOINT=<your-endpoint-uuid>
python scripts/download_forcing.py --ocean        # thetao/so/tf + climatology + bias
python scripts/download_forcing.py --calibration  # meltMIP obs melt, IMBIE2 basins, grid, topography
python scripts/download_forcing.py --scenarios    # per-(ESM, scenario) forcing (cores 1-8)
python scripts/download_forcing.py --scenarios --esm MRI-ESM2-0 --scenario historical,ssp585
python scripts/download_forcing.py --status
```

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_GLOBUS_COLLECTION` | source collection UUID | `ccc9bbd2-4091-4e35-addd-eeb639cf5332` |
| `GLOBUS_LOCAL_ENDPOINT` | your Globus Connect Personal endpoint UUID | required to transfer |

Scenario forcing lives in the collection's top-level `/ISMIP7/AIS/<ESM>/<scenario>/`
tree. `--scenarios` mirrors the minimal runtime sets (SDBN1-8000m `acabf` and
`acabf-anomaly`, ocean `tf` and `so`, fracture) with version autodetection and
checksum sync, so re-runs are completeness checks. Climatology, obs and
calibration sets come from `/ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/AIS`; the
`OCEAN_FILES` and `CALIBRATION_FILES` dicts in the script are the manifest.
Without `GLOBUS_LOCAL_ENDPOINT` the script prints the paths for a manual
transfer in the web app.

### 2b. The runtime tree

```
ISMIP7/AIS/
  <ESM>/<scenario>/<SDBN1|GEMB-SDBN1>-8000m/<var>/<version>/   # acabf, acabf-anomaly, ts, tas, pr
      <var>_AIS_<ESM>_<scenario>_<product>_<version>_<YEAR>.nc
  <ESM>/<scenario>/ocean/<tf|thetao|so>/<version>/
  <ESM>/<scenario>/fracture/[v*/]                     # collapse and lake masks
  meltMIP/OI_Climatology_ismip8km_60m_<tf|so|thetao>_extrap.nc
  parameterisations/ocean/imbie2/                     # basin numbers for per-basin K
  parameterisations/{ocean,fracture}/
```

Readers pin `version=v2` for the atmosphere and `v3` for the ocean, falling
back to the highest `v<N>` present, dotted versions included, so the `v2.1`
fracture release and MRI-ESM2-0's `v1` resolve without code changes. The
atmosphere directory is whichever of `SDBN1-8000m` and `GEMB-SDBN1-8000m`
exists, `SDBN1` first, so trees fetched before MRI's August 2026 rename still
work. Fracture masks resolve flat or versioned.

**One year bridges the end of a series.** CESM2-WACCM's atmosphere stops at
2299 and the empty 2300 files were withdrawn, so a request for the single year
after the last one on disk reuses that year and logs it once per variable.
Anything further past the end, or a gap inside the series, raises.

**Per-year atmosphere files hold 12 monthly slices.** The reader collapses them
to the annual mean, weighting months by length from `time_bnds`, falling back
to coordinate spacing and then to an unweighted mean with a warning. A
length-1 time axis passes through unchanged.

API: `ISMIP7Atmosphere(esm, scenario).get_smb(year, x, y, anomaly=...)`,
`ISMIP7Ocean(...).get_thermal_forcing(...)` and `.get_salinity(...)`,
`ISMIP7Fracture(...).get_collapse_mask(year, x, y)`, and
`make_forcing_callback(atm=, ocean=, fracture=, K=, K_per_basin_npz=)` which
bundles them into the callback `run_simulation` expects.

---

## 3. Build the mesh

Adaptive isotropic mesh with Ua-style grounding-zone and calving-front
refinement, sized from BedMachine geometry and MEaSUREs strain rate.

```bash
cd antarctica
python scripts/mesh_antarctica.py --lc 2500 --lc-coarse 64000 --buffer-m 20000
ISMIP7_LC=2500 ISMIP7_LC_COARSE=64000 ISMIP7_BUFFER_M=20000 python scripts/mesh_antarctica.py
# dev mesh for inversion_icepack2.py, diagnostic_solve.py and run_eigendec.py:
python scripts/mesh_antarctica.py --lc 8000 --lc-coarse 80000 --buffer-m 20000
# mesh/antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.msh
# mesh/boundary_ids_antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.json
```

`--lc` and `--lc-coarse` are the fine (grounding line, calving front) and
coarse (interior) element sizes in metres. GL-band sizes, `shelf_size`,
`buffer_size` and the strain-rate floor scale with `lc/2500`; calving-front
decay lengths are floored at `lc` and `1.25*lc`.

`--buffer-m` pushes the outline into the ocean before meshing (default 20000,
`0` for none), letting icepack2 handle `h=0` at an interior calving front
instead of a `calving_terminus` BC. Both the `.msh` and its sidecar are named
after the exact `(COARSE, FINE, BUFFER_M)` triple, so builds never collide.

**Boundary ids.** gmsh physical groups alternate `Calving_0, Other_1, ...`,
auto-numbered from 1, so odd tags are calving and even tags are other.
`mesh_antarctica.py` writes that split to the per-mesh sidecar (named by
`scripts/mesh_naming.py`). To regenerate one without rebuilding the mesh:
`ISMIP7_BUFFER_M=<N> python scripts/make_boundary_ids.py`.

Solvers resolve the sidecar through `icepack2_tools/boundary.py`:
`ISMIP7_BNDIDS`, then the per-mesh `mesh/boundary_ids_<mesh stem>.json`, then
the shared `mesh/boundary_ids.json`. Readers hard-error on an unclassified
exterior marker or a missing id. A silent mismatch once left 95% of the ice
front without calving back-pressure. Runs print the covered front length in km
and percent.

A forward takes the sidecar name from the MAP or restart checkpoint, which
records `mesh_basename` and the `lc`, `lc_coarse` and `buffer_m` parameters, so
environment drift cannot swap a sidecar mid-trajectory.

---

## 4. Invert for basal and rheology fields (`inversion_icepack2.py`)

MAP estimate of bed friction `θ` and rheology `φ` from the diagnostic 3-field
(V x Σ x τ) system, regularized, with n=1 to 3 continuation, through
`tlm_adjoint`.

```bash
cd antarctica
ISMIP7_LC=2500 mpiexec -n 12 python scripts/inversion_icepack2.py
# mesh/inversion_icepack2_<budd|rc>_n3_dg0_<LC>.h5
```

The controls are log deviations from physical priors: `θ = log(C/C_w0)` on the
balance-friction anchor and `φ = log(A/A_prior)` on a thermomechanical fluidity
prior computed at setup and stored in the MAP. That prior reads the section 1
RACMO SMB and the section 2 `tas` climatology, so download the forcing first.
See `N3_FRAMEWORK.md`, and `ISMIP7_FLUIDITY_PRIOR=legacy` to skip it.

`ISMIP7_FRICTION` selects the law (`budd` or `regularized_coulomb`, tagged
`_budd` and `_rc` in the filename). The MAP name carries the law, the flow
exponent and the geometry space, built by `icepack2_tools/naming.py`, which the
forward, `preflight.py` and the gates all import. A MAP is valid only for its
own geometry space, since the inversion absorbs the calving-front treatment into
`θ` and `φ`. A forward that finds only a legacy untagged MAP loads it with a
loud warning. See
`../GEOMETRY_DISCRETIZATION.md`.

### Transient (dH/dt-constrained) inversion

A velocity-only inversion fits `u` while leaving `div(h u)` unconstrained, so
the MAP can carry a flux divergence inconsistent with the observed geometry.
`ISMIP7_DHDT_WEIGHT > 0` adds one implicit-Euler prognostic step after the
diagnostic solve, using the model's own DG0 upwind operator, and scores the
resulting tendency against the observed mean dH/dt
(`icepack2_tools/obs_dhdt.py`, from the MIPkit). It needs
`ISMIP7_GEOMETRY_SPACE=dg0` and applies to grounded ice.

```bash
ISMIP7_DHDT_WEIGHT=1.0 ISMIP7_MAP_OUT=mesh/inversion_transient_2500.h5 \
  ISMIP7_LC=2500 mpiexec -n 12 python scripts/inversion_icepack2.py
python scripts/compare_dhdt.py vel=mesh/<velocity-only>.h5 tr=mesh/<transient>.h5
```

`compare_dhdt.py` is the payoff diagnostic: both MAPs fit `u`, so the thickness
tendency is the observable that separates them. It drives the step with the
`velocity` stored in the MAP. The inversion saves that field only when the
final solve converged at the MAP's own controls and its misfit agrees with the
last accepted optimization state. MAPs written before that guard can carry a
bad one, and a score built on it means nothing. The controls are unaffected:
a forward re-solves the diagnostic from `θ` and `φ`.

`ISMIP7_DHDT_NET_SIGMA > 0` adds a term on the integrated grounded dH/dt. It is
off by default so the integrated trend stays independent validation against
IMBIE and GRACE. The per-iteration `net=` diagnostic prints either way. See
`reports/ISSUE_DRAFT_net_mass_balance_term.md`.

The MAP filename encodes friction, `LC`, geometry space and flow exponent only,
so velocity-only and transient variants collide. Give variants their own
`ISMIP7_MAP_OUT`. Every MAP records its objective as root attributes
(`misfit_norm`, `gamma_theta`, `gamma_phi`, `log_vel_weight`, `log_vel_eps`,
`dhdt_weight`, `dhdt_net_sigma`, `mesh_basename`, `lc`, `lc_coarse`,
`buffer_m`):

```bash
python -c "import h5py,sys; print(dict(h5py.File(sys.argv[1])['/'].attrs))" MAP.h5
```

---

## 5. Calibrate ocean melt (`calibrate_melt.py`)

Solves for the Burgard quadratic-mixed-slope coefficient K, global and per
IMBIE2 basin, against integrated observed shelf melt. The target is the July
2026 table combining Paolo, Davison and Adusumilli, 1067.4 Gt/yr, read from
`<DATA_ROOT>/meltobs/Melt_Paolo_Davison_Adusumilli_imbie2.csv`. With that file
absent it falls back to the older Paolo and Adusumilli table (865.0 Gt/yr) under
`<DATA_ROOT>/parameterisations/ocean/meltobs/`; `ISMIP7_MELT_OBS_CSV` names
either. Needs section 2 forcing and a section 4 mesh.

The newer table comes from the Source Cooperative melt-calibration product:

```bash
python antarctica/scripts/download_mirror.py \
    --product ismip7-ais-melt-calibration data/meltobs/
```

```bash
cd antarctica
ISMIP7_LC=2500 python scripts/calibrate_melt.py
# results/calibrated_K_per_basin_<LC>.npz
```

Only the mesh is read from the MAP, so any MAP built on it serves. The default
is the section 4 name for the configured `ISMIP7_FRICTION`; `ISMIP7_INV_H5`
names a different one.

The control requires this npz. Projections take it (`K_per_basin_npz=`) or a
scalar `ISMIP7_K_MELT`.

---

## 6. Control and projections

Forward runs go through `scripts/simulation.py` (`setup_model` and
`run_simulation`). Drivers live in `scripts/control/`, `scripts/historical/`
and `scripts/projections/`.

```bash
cd antarctica
mpiexec -n 12 python scripts/control/run.py
# results/ctrl2015_<esm>_<lc>_{final.h5, t<year>.h5, timeseries.csv}
mpiexec -n 12 python scripts/projections/ssp585_cesm_waccm.py
# results/ssp585_cesm2_waccm_<lc>_{final.h5, timeseries.csv}
```

The other scenario drivers (ssp126 and ssp370 for both ESMs, plus `ocx.py`) are
shims over `scripts/experiment.py`. Run a historical driver first to produce
`results/hist_<esm>_<lc>_final.h5`: projections and the control both branch
from it, so they share a t=0 state and the same frozen apparent-MB correction,
and their relaxation drift cancels in projection minus control (the ISMIP6
ctrl_proj convention). Without it a projection cold-starts from BedMachine and
the control warns that it starts from a different geometry.

Run management on the drivers: `--restart <ckpt>` or `ISMIP7_RESTART` resumes;
`ISMIP7_AUTO_RESUME=1` picks up the newest checkpoint for the experiment, which
is what lets a chained batch job continue itself, and takes precedence over the
historical endpoint so only the first link starts there; `--tag` or
`ISMIP7_RUN_TAG` suffixes the experiment name so a method line keeps and
resumes its own files; `--checkpoint-interval` sets the step-count fallback.
Checkpoints carry the mesh, geometry, inversion fields and the full `(u, M, τ)`
state, so restarts work at any rank count. A resume refuses to start when
`ISMIP7_FRICTION` or `ISMIP7_APPARENT_MB` disagree with the checkpoint.

**Is the run on track?**

```bash
python scripts/check_ismip6_track.py results/<exp>_timeseries.csv
# exit code 0 when no FAIL rows, so gates can chain on it
```

It audits a timeseries against IMBIE dM/dt, Rignot melt and calving, RACMO SMB
and ISMIP6-class control drift, with a runaway detector.
`scripts/compare_ismip6.py <proj.csv> <ctrl.csv>` overlays projection minus
control sea-level contribution on the ISMIP6 ensemble.

Read-only diagnostics:

| | |
|---|---|
| `compare_runs.py LABEL=results/<a> ...` | overlays budget timeseries to show where two runs part ways |
| `plot_movie.py results/<exp>` | thickness change, speed and thickness frames plus an mp4 under `figs/movie_<exp>/` |
| `region_budget.py <ckpt>.h5 [<later>.h5] [--csv <run>_timeseries.csv]` | splits the budget into grounded and floating ice, so a control that gains volume above flotation is read against the observed 2000 to 2200 Gt/yr of discharge |
| `score_map.py MAP.h5 [...]` | scores an inversion by `Q(u_model)/Q(u_obs)` across its own grounding line, overall and per speed band; re-solves through `setup_model`, so periodic MAPs without a velocity work |
| `plot_map.py MAP.h5 [--diff B.h5]` | model and observed speed and their difference, `θ`, `C = C_w0 exp(θ)` on grounded ice, and `φ`, into `figs/maps/` |

`region_budget.py` and `score_map.py` take the run's environment, which must
match the inversion's.

### Calving front on a buffered mesh (`ISMIP7_CALVING`)

A buffered mesh has no calving sink: ice reaching the 2015 outline flows into
empty buffer cells and only shelf melt removes mass, so a control gains mass
(+1000 to +1500 Gt/yr in the September 2026 32 km case with
`ISMIP7_APPARENT_MB` unset). `ISMIP7_FIXED_FRONT` removes what crosses the
outline while leaving it fixed.

`ISMIP7_CALVING` replaces that with the shared level set
(`icepack_tools.levelset`, wrapped by `icepack2_tools/levelset.py`, the object
CalvingMIP also runs). Each step `phi` solves the eikonal problem
`|grad phi| = 1` with `phi = 0` on the facets between ice and ice-free cells of
the transport's own thickness, negative in ice and positive in water. The
calving rate `c` then retreats the front by `phi_t - c|grad phi| = 0` (Hahn,
Mikula and Frolkovic 2025, arXiv:2504.05845), linearised with the previous unit
gradient and solved by cell-centred finite volumes.

Advance needs no extra mechanism: the upwind DG0 transport fills any cell the
ice flows into and the next extent includes it. Removal conserves calved mass
in three parts, all booked to the `calv` column:

1. cells the front passed entirely (`phi > 0`) are emptied;
2. every front cell sheds `min(1, c dt L/A)` of its thickness, the mass
   `c h L dt` a front retreating at `c` loses, which carries sub-cell retreat
   between steps;
3. the sliver left in a cell that held ice when the step began.

A cell below `ISMIP7_FRONT_HMIN` (1 m) that already held ice is a retreating
front cell and is emptied. A cell that was ice-free keeps whatever the
transport put there, which is how the front advances. Outside the extent the
thickness can therefore be small and nonzero.

Momentum needs no front term: under DG0 geometry the facet term
`rho g avg(h) jump(s)` at an ice/water face is the terminus water-pressure
force. The one momentum change is the drag gate: the mask is 1 only where `phi`
exceeds one cell diameter, so floor-cell ocean drag acts only in water further
than a cell from the front. Thin cells inside the t=0 extent are damped by
`h_visc_floor`, the friction law and the `ISMIP7_ALPHA_GL` collar. In the strip
a free law advances into, `C_w0` comes from the t=0 geometry where `H = 0`, so
`tau_b = 0` under both laws and the damping is `h_visc_floor` and the collar.

`vonmises` is Morlighem et al. 2016 verbatim:
`c = |u| sqrt(3) B eps~^(1/n) / sigma_max`, with `eps~` from the tensile
principal strain rates, `B = A^(-1/n)`, and separate grounded and floating
thresholds. Those thresholds are the tuning targets: a 2015 control should hold
the observed front (the obs kit's 24 yearly Greene masks, 1997 to 2021) and
discharge about 1300 Gt/yr. The level set is checkpointed as `levelset` for
diagnostics; a restart rebuilds the front from the thickness. The exception is
`fixed`, which anchors on `H_init` so a resumed run does not re-freeze the
front where it restarted. The shared implementation's tests are
`icepack_tools/test/levelset_test.py`; the ISMIP7-side rules (retreat-sliver
mask, apparent-MB extent masking, the `fixed` law's t=0 anchor) are covered by
`tests/`. The level-set unit tests written against this integration in Sep 2026
were lost before they were committed and are still to be rebuilt.

**Control and projection configurations differ.** The protocol's control is an
unforced constant-climate run with calving set to end-of-2014 conditions, so
the control here is `ISMIP7_APPARENT_MB` with `ISMIP7_FIXED_FRONT=1` and
`ISMIP7_CALVING=none`. That is what `run_core_matrix.sh` runs and what every
control result used. `ISMIP7_CALVING=fixed` also pins the front and is a
different run: it builds a level set, so ocean drag is gated off near the front
and the retreat-sliver rule applies inside the t=0 extent, giving a slightly
different `calv` column and settled front. Both close the budget.

`vonmises` is for projections. A configured law owns the front outright, so the
legacy `ISMIP7_FIXED_FRONT` mask removes nothing when a law is set and
`vonmises` is never silently pinned. The apparent-MB reference `a_ref` is
defined only on the t=0 ice extent under every law. Under a free law it is also
cleared each step wherever the level set reports ice-free, irreversibly, so a
calved cell is not regrown and an advanced-into cell is not re-emptied. The run
log prints one `Calving front owner:` line naming the mechanism in force.

### The whole matrix in one command (`run_core_matrix.sh`)

Runs cores 1 to 11 in dependency order (historicals, controls, projections,
OCX) and finishes with the audit and ensemble comparison for each core that
reached its target year.

```bash
ISMIP7_RUN_TAG=n3 ISMIP7_LC=32000 antarctica/scripts/run_core_matrix.sh
CORES=1,2,9 antarctica/scripts/run_core_matrix.sh
```

Cores run sequentially, load-gated by `MAX_LOAD` (default cores minus 8). A
core already at its target year is skipped. Output predating the annual-mean
forcing fix is archived under `results/archive_stale_<stamp>/` and re-run.
The runner's own knobs (`CORES`, `MAX_LOAD`, `MAX_ATTEMPTS`, `FRESH`, `REUSE`,
`NRANKS`, `PROV_REF`) and the reuse rules are documented in its header.

Record each completed core with `python scripts/core_report.py --core <N>
--name <exp> --csv <timeseries.csv> --log <run.log>` (add `--ctrl-csv` for a
projection), run in that run's own shell so it captures the environment.
`--superseded "<reason>"` stamps a record when a later run replaces it.

### Environment knobs (inversion)

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_MAP_OUT` | output path for the MAP, overriding the generated name. Use it for smoke tests and variants so a short run cannot replace a production MAP. A bare filename resolves under `mesh/` | generated |
| `ISMIP7_MISFIT_NORM` | `sigma` divides each residual by its datum's squared error, giving a dimensionless chi^2; `none` is the legacy dimensional misfit. Selects the `ISMIP7_GAMMA_*` defaults | `sigma` |
| `ISMIP7_LOG_VEL_WEIGHT` | weight on the ISSM logarithmic velocity misfit (cost function 103). The chi^2 alone over-weights slow interior ice and leaves discharge-carrying tributaries 40 to 50% too slow; the log term is scale free. `auto` equalises it with the chi^2 term at the warm-start state. Stamped into the MAP | `0` |
| `ISMIP7_LOG_VEL_EPS` | regularisation speed (m/yr) inside the log | `1.0` |
| `ISMIP7_GAMMA_THETA` / `ISMIP7_GAMMA_PHI` | Whittle-Matern prior strength on `θ` and `φ`, coupled to `ISMIP7_MISFIT_NORM` since normalising divides the misfit by about sigma^2 | `1e5` under `sigma`, `1e4` under `none` |
| `ISMIP7_L_REG` | prior correlation length (m) | `7.5e3` |
| `ISMIP7_MAXITER` | L-BFGS-B iteration cap | `500` |
| `ISMIP7_GRAD_PRECOND` | `none` is the raw-dof l2 metric, which is mesh dependent, so fine grounding-line cells converge slowest. `mass` optimises in `u = sqrt(M) x`, making the rate mesh independent. Defaults to `none` to keep runs comparable with everything measured so far | `none` |
| `ISMIP7_SIGMA_U_FLOOR` | floor on the per-component MEaSUREs error (m/yr), so near-zero errors cannot let a few nodes dominate | `1.0` |
| `ISMIP7_SIGMA_U_UNOBS` | sigma (m/yr) where MEaSUREs reports no error. Those nodes carry a zero-filled `u_obs`, so they need a large sigma when `ISMIP7_OBS_MASK=0` | `1e4` |
| `ISMIP7_OBS_MASK` | `0` drops the velocity-observation mask | `1` |
| `ISMIP7_DHDT_WEIGHT` | weight on the dH/dt chi^2; `0` disables the transient constraint. Needs `ISMIP7_GEOMETRY_SPACE=dg0` | `0` |
| `ISMIP7_DHDT_SIGMA` | assumed dH/dt uncertainty (m/yr), hand set since the MIPkit ships no uncertainty field | `0.1` |
| `ISMIP7_DHDT_DT` | timestep of the prognostic step (yr) | `1.0` |
| `ISMIP7_DHDT_VAR` | `dhdt_smith` (firn corrected, 2003 to 2019) or `dhdt_cpom` (uncorrected, so not interchangeable) | `dhdt_smith` |
| `ISMIP7_DHDT_MELT` | `0` drops ocean melt from the prognostic source. Per-basin K comes from `ISMIP7_K_PER_BASIN_NPZ`, else the `<lc>` npz, else the 2500 m file; with none it warns and uses SMB only, melt being zero on grounded ice | `1` |
| `ISMIP7_DHDT_CLIM_START` / `_END` | RACMO climatology window for that source | `2003` / `2019` |
| `ISMIP7_DHDT_REACH` | pixel-to-cell reach as a multiple of `sqrt(area)`, rejecting pixels outside the mesh that nearest-centroid assignment would snap onto boundary cells | `0.75` |
| `ISMIP7_DHDT_NET_SIGMA` | sigma (Gt/yr) on the integrated grounded dH/dt; `0` disables the net term. Active only with `ISMIP7_DHDT_WEIGHT > 0` | `0` |
| `ISMIP7_OBS_KIT` | path to `AntarcticaObsISMIP7-v*.nc`. The kit is needed only to build the dH/dt cache rasters in `antarctica/data/dhdt_cache/`; with those staged it may be absent. A path that does not exist is a hard error | newest under `<DATA_ROOT>/obs/mipkit` |

### Environment knobs (forward runs)

`icepack2_tools/runconfig.py` owns the run-shaping knobs (`ISMIP7_LC`,
`ISMIP7_LC_COARSE`, `ISMIP7_FRICTION`, `ISMIP7_GEOMETRY_SPACE`,
`ISMIP7_N_FLOW`), so an unset knob cannot mean one resolution to the inversion
and another to the preflight.

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_LC` / `ISMIP7_LC_COARSE` | fine and coarse mesh resolution tags, selecting mesh and MAP | `2500` / `64000` |
| `ISMIP7_BUFFER_M` | outline buffer (m) in the default mesh and sidecar names | `20000` |
| `ISMIP7_MESH` | mesh path for the inversion and tools. A forward takes its mesh from the checkpoint | `mesh/antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.msh` |
| `ISMIP7_RASTER_SAMPLE` | how BedMachine lands on a DG0 cell. `vertex` projects the CG1 vertex interpolant; `cell_mean` takes the raster's true cell mean. `cell_mean` measured rougher: neighbouring cells share two of three vertex samples, so `vertex` damps jumps by construction. Cell means raised interior surface jumps 6% and bed and thickness jumps 35%, and at 2 km the momentum solve did not converge within 60 minutes. It does classify flotation better (32 km misclassification 9.1% to 3.2%), so the knob stays. Stamped into the MAP and read back by the forward. Reproduce with `probe_raster_sampling.py` | `vertex` |
| `ISMIP7_INVERSION` | explicit MAP path for a forward or preflight. The forward checks the MAP's recorded `friction`, `n_flow` and `geometry_space` against the run and aborts on a mismatch, warning only when the MAP predates those attributes; `preflight.py` checks that the file exists. Use it to A/B MAPs on one mesh | derived |
| `ISMIP7_CALVING` | `none`, `fixed` or `vonmises` (see above) | `none` |
| `ISMIP7_CALVING_SIGMA_MAX_GROUNDED` / `_FLOATING` | von Mises thresholds (MPa) | `1.0` / `0.15` |
| `ISMIP7_FRACTURE` | `mask` applies the ISMIP7 collapse forcing to floating cells, booked as calving. Masks exist for the SSPs only, so the control, historicals and OCX abort on `mask`. Needs DG0 | `none` |
| `ISMIP7_OUTPUT` | `1` records the ISMIP7 yearly fields and scalars (`<exp>_<lc>_ismip7_annual_<year>.h5`, `<exp>_<lc>_ismip7_scalars.csv`), regridded afterwards by `write_ismip7_output.py`. The value set is closed, so a typo is rejected at startup. A chained projection must export it on every link; a link that cold-starts mid-year logs the gap and begins at the next 1 January. Resuming continues a series, and a cold start into a populated series is refused | unset |
| `ISMIP7_BNDIDS` | boundary-id JSON override | per-mesh sidecar, else `mesh/boundary_ids.json` |
| `ISMIP7_GEOMETRY_SPACE` | `dg0` (one thickness for terminus force and mass flux) or `cg1` (legacy, A/B only). Selects the MAP. See `../GEOMETRY_DISCRETIZATION.md` | `dg0` |
| `ISMIP7_DATA_ROOT` | forcing tree root | `<repo>/ISMIP7/AIS` |
| `ISMIP7_T_END` / `ISMIP7_DT` | end time and timestep (yr). `t=Y.0` is 1 January of year Y, so a run covering 2015 to 2300 ends at `2301` and a historical covering 1850 to 2014 ends at `2015`. Each driver owns its end (historical `2015`, ssp370 `2101`, other projections and control `2301`, OCX `2026`) | driver's own / `1.0` |
| `ISMIP7_FRICTION` | `budd` or `regularized_coulomb`; selects the MAP | `budd` |
| `ISMIP7_OUTPUT_INTERVAL` | timeseries row every N steps | `10` |
| `ISMIP7_CHECKPOINT_EVERY_YR` / `ISMIP7_KEEP_CHECKPOINTS` | checkpoint cadence in model years, and how many to keep besides `_final.h5` | `5` / `3` |
| `ISMIP7_RESTART` | restart checkpoint | `hist_<esm>[_<tag>]_<lc>_final.h5` if present |
| `ISMIP7_AUTO_RESUME` | resume from this experiment's newest checkpoint when no `ISMIP7_RESTART` is given. An integer flag, `=0` disables it, since the runners export it unconditionally and `--export=ALL` cannot unset. `projection.sbatch` refuses to chain when it is off | unset |
| `ISMIP7_RUN_TAG` | experiment-name suffix for a parallel method line | unset |
| `ISMIP7_WALL_STOP_MIN` | wall-clock budget in minutes from process start, checked before each step against the longest step so far, so the run writes its final checkpoint and exits with `t_yr` short of `t_end` for a chained job to resume. `projection.sbatch` derives it from the job's own TimeLimit, holding back 25 minutes. `0` disables | `0` |
| `ISMIP7_EXPERIMENT_NAME` | the run's identity, used by `adapt_mesh.py` to name adapted meshes and sidecars so parallel experiments cannot overwrite each other. Set by `run_adaptive.py --experiment-name`. See `../UA_ADAPTIVE_MESH.md` | unset |
| `ISMIP7_APPARENT_MB` | `1` or `balance` zeroes the t=0 thickness tendency (ISMIP6 ctrl_proj style); `div` cancels only the flux divergence; `0`, `off`, `none` and empty disable it | unset |
| `ISMIP7_FIXED_FRONT` | hold the calving front at the t=0 extent, tallying inflow beyond it as calving. `=0` disables. Ignored whenever an `ISMIP7_CALVING` law is configured | unset |
| `ISMIP7_LEGACY_TRANSPORT` | restore the pre-July-2026 CG-projection transport (needs `cg1`) | unset |
| `ISMIP7_SNES_TYPE` / `ISMIP7_SNES_MAXIT` | diagnostic Newton type and iteration cap | `newtonls` / `200` |
| `ISMIP7_K_MELT` / `ISMIP7_K_PER_BASIN_NPZ` | scalar Burgard K (projections), per-basin K file (control) | `1.15e-4` / `results/calibrated_K_per_basin_<lc>.npz` |
| `ISMIP7_MELT_OBS_CSV` | per-basin melt observation table read by `scripts/calibrate_melt.py`; columns are located by header name, so either published table serves | `<DATA_ROOT>/meltobs/Melt_Paolo_Davison_Adusumilli_imbie2.csv`, else the older Paolo and Adusumilli table |
| `ISMIP7_ESM` | ESM for the control | `CESM2-WACCM` |
| `ISMIP7_CLIM_SCENARIO` / `_START` / `_END` | reference-climate pool: the scenario pooled with `historical`, and the window, shared by the control's SMB climatology and the projections' aSMB re-reference through `icepack2_tools/climatology.py`. A partial pool warns | `ssp126` / `2000` / `2029` |
| `ISMIP7_H_CLAMP` | thickness floor (m) | `0` |
| `ISMIP7_NO_CALVING_TERMINUS` | drop the calving-terminus BC | unset |
| `ISMIP7_SUBCYCLES` / `ISMIP7_RESCUE_MAXIT` | dt-subcycle rescue ladder, and the Newton cap on its rungs | `1,4,16` / `600` |
| `ISMIP7_H_OCEAN` / `ISMIP7_K_LIM` | front backstops read by `scripts/simulation.py`: the thickness (m) at which the ice-free ocean drag ramps to zero, and the speed-limiter coefficient the rescue ladder raises for a rescue solve | `10.0` / `1e-3` |
| `ISMIP7_ALPHA_GL` | grounding-line coercivity, Budd only, read by `scripts/simulation.py` and `scripts/inversion_icepack2.py` | `0.5` (`0` for RC) |
| `ISMIP7_RC_HVISC_FLOOR` / `ISMIP7_RC_CW0_FLOOR` | RC viscous-thickness and `C_w0` floors, read by `scripts/simulation.py` and `scripts/inversion_icepack2.py` | `10.0` / `0.0` |
| `ISMIP7_M_SLIDE` | sliding exponent, read by `scripts/inversion_icepack2.py`, `scripts/simulation.py`, `scripts/thermo_prior.py`, `scripts/plot_map.py`, `scripts/run_eigendec.py` | `3.0` |

> **dt guidance.** Use `ISMIP7_DT=0.1` for production projections. `0.25` is
> acceptable when 10 steps per year is too costly. `dt=1.0` mis-melts per step
> and resurrects clamped cells.

---

## Outputs

Per experiment in `results/`:

- `<exp>_final.h5`, a self-contained restart checkpoint: mesh, geometry,
  inversion fields, the full `(u, M, τ)` state, the frozen apparent-MB
  reference, and `levelset` when a calving law is configured. Under `dg0` the
  saved `thickness` is the prognostic state; a `cg1` run also saves
  `thickness_dg`. The `geometry_space` and `mesh_basename` attributes let a
  restart resolve the same sidecar.
- `<exp>_t<year>.h5`, periodic checkpoints, keeping the
  `ISMIP7_KEEP_CHECKPOINTS` most recently written.
- `<exp>_timeseries.csv`, one row per `OUTPUT_INTERVAL` steps:
  `year, vaf_mm_sle, mass_gt, smb_gtyr, melt_gtyr, outflux_gtyr, calv_gt,
  clamp_gt, resid_gt, amb_gtyr`. The residual must close to 0.00.

VAF is in mm of sea-level equivalent, mass in Gt.

---

## Known issues

**Diagnostic-Newton wall on hard projection geometries.** The forward blow-ups
are fixed (balanced apparent-MB init plus persistent DG0 thickness state), and
the Newton can still stall on evolved geometries. `ISMIP7_SNES_TYPE` and
`ISMIP7_SNES_MAXIT` are the knobs. The aSMB-forced walls seen so far predate
the annual-mean forcing fix and may be forcing induced; see
`reports/MATRIX_STATUS.md`.

**Wall retry.** When the rescue ladder is exhausted the run saves and stops
short of its target year. Relaunching from that state has cleared the wall in 3
of 3 observed cases (ssp585-CESM at 2096.7, CTRL-CESM at 2268, CTRL-MRI at
2250, both CTRLs then reaching 2300), because a fresh process re-runs the n=1
to n continuation at the loaded geometry. `run_core_matrix.sh` does this
automatically on a workstation. On a cluster the relaunch is yours: resubmit,
and `ISMIP7_AUTO_RESUME` picks the run up. A checkpoint written with
`stalled=1` makes `projection.sbatch` report and exit 1 without a successor, so
an unattended chain cannot spend days re-attempting the same years.

**Mesh and sidecar naming is per `(COARSE, FINE, BUFFER_M)`.** Older meshes
built before this convention need renaming or rebuilding.

**`icepack2_tools/coupled.py`** sketches ice-plume coupling against a
`PlumeModel` that does not exist in this tree. Nothing runs it.

---

## References

- Burgard et al. 2022, *The Cryosphere*, basal-melt parameterisation assessment.
- multimelt reference implementation: https://github.com/ClimateClara/multimelt
- ISMIP7 ocean forcing pipeline: https://github.com/ismip/ismip7-antarctic-ocean-forcing
- Greenland companion: https://github.com/dlilien/ISMIP7_Greenland_Icepack
