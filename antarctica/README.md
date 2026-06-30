# ISMIP7 Antarctica — data, setup, and how to run

Antarctic ISMIP7 submission built on **icepack2 / Firedrake**. This document is
the end-to-end guide: what to install, how to download every input (including the
Globus forcing), and the exact command sequence to reproduce a run. Companion to
David Lilien's Greenland repo (https://github.com/dlilien/ISMIP7_Greenland_Icepack).

> New here? Read top to bottom once. The pipeline is a strict dependency chain —
> each step consumes the previous step's output (data → mesh → inversion →
> calibration → control/projection).

---

## 0. Prerequisites

### Software

- **Firedrake** with **icepack2** (the mixed/3-field formulation) and
  `tlm_adjoint` for the inversions. Firedrake brings its own PETSc/petsc4py/MPI.
  Activate that environment before anything else, e.g.:
  ```bash
  source ~/venv-firedrake/bin/activate
  ```
- Python packages used by the scripts (install into the Firedrake venv):
  ```bash
  pip install globus-sdk earthaccess xarray netCDF4 scipy shapely gmsh
  ```
- **gmsh** (the `gmsh` Python module above is sufficient for meshing).

### Accounts / credentials

| For | Account | Where |
|-----|---------|-------|
| BedMachine, MEaSUREs velocity (NSIDC) | NASA Earthdata (free) | https://urs.earthdata.nasa.gov/users/new |
| ISMIP7 ocean/atmosphere forcing | Globus + access to the ISMIP6/7 collection | https://app.globus.org |

To pull data **to this machine** over Globus you also need **Globus Connect
Personal** running locally and its endpoint UUID
(https://www.globus.org/globus-connect-personal).

### Repo layout (what lives where)

```
antarctica/
  data/            # observational inputs (BedMachine, velocity, RACMO)  [gitignored]
  mesh/            # *.msh + boundary_ids_antarctica_*.json + inversion_*.h5  [gitignored except boundary_ids*.json]
  results/         # checkpoints (*.h5), timeseries (*.csv), logs        [gitignored]
  scripts/         # all entry points (see §3–§6)
ISMIP7/AIS/        # ISMIP7 forcing tree the *runtime* reads             [gitignored]
icepack2_tools/    # reusable library (mesh, forcing, eikonal, grounding, regrid)
```

Everything large is gitignored. The only tracked files under `antarctica/mesh/`
are the tiny `boundary_ids_antarctica_*.json` (the gmsh physical-line →
calving/other map; see §3).

---

## 1. Observational data (`download_data.py`)

BedMachine geometry, MEaSUREs velocity, and RACMO SMB. RACMO is public (Zenodo);
the NSIDC products need an Earthdata login (handled interactively by
`earthaccess`).

```bash
cd antarctica
python scripts/download_data.py
```

| Dataset | Product | Auth | → lands in |
|---------|---------|------|-----------|
| BedMachine Antarctica v4 | NSIDC-0756 | Earthdata | `data/bedmachine/` |
| MEaSUREs Ice Velocity v2 | NSIDC-0484 | Earthdata | `data/velocity/` |
| RACMO2.4p1 SMB | Zenodo `10.5281/zenodo.14217231` | none | `data/racmo/` |

The scripts skip files that already exist, so re-running is cheap.

---

## 2. ISMIP7 forcing via Globus

There are **two distinct things** here — read this carefully, it's the most
common source of confusion:

1. **The runtime tree `ISMIP7/AIS/`** — this is what the simulation actually
   reads at run time (via `icepack2_tools/forcing.py`). Its layout follows the
   official ISMIP7 protocol (see §2b). Point the code at it with
   `ISMIP7_DATA_ROOT` (defaults to `<repo>/ISMIP7/AIS`).
2. **`download_forcing.py`** — a helper that pulls the *ocean* CMIP/climatology/
   calibration files from the ISMIP6/7 Globus collection's
   `share_with_modellers/` area into `data/forcing/`. Its directory layout is the
   `share_with_modellers` one, **not** the `ISMIP7/AIS/` runtime layout — so files
   fetched this way currently need to be reorganized/symlinked into `ISMIP7/AIS/`
   (or pointed at via `ISMIP7_DATA_ROOT`). See the caveat at the end of §2a.

### 2a. Using `download_forcing.py`

```bash
cd antarctica

# one-time: authenticate (opens a Globus URL, paste back the auth code).
# Token is cached at ~/.ismip7_globus_tokens.json (chmod 600).
python scripts/download_forcing.py --login

# browse the remote tree to sanity-check paths
python scripts/download_forcing.py --list

# tell the script where to *put* the files (your Globus Connect Personal endpoint)
export GLOBUS_LOCAL_ENDPOINT=<your-endpoint-uuid>

# download (omit a flag to get everything; --dry-run to preview)
python scripts/download_forcing.py --ocean        # CESM2-WACCM thetao/so/tf + climatology + bias
python scripts/download_forcing.py --calibration  # meltMIP obs melt, IMBIE2 basins, grid, topography
python scripts/download_forcing.py --status        # what's present locally
```

Knobs:

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_GLOBUS_COLLECTION` | source collection UUID | `ccc9bbd2-4091-4e35-addd-eeb639cf5332` |
| `GLOBUS_LOCAL_ENDPOINT` | **your** local Globus Connect Personal endpoint UUID | _(required to transfer)_ |

Remote base path on the collection:
`/ISMIP6/ISMIP7_Prep/AIS_ocean/share_with_modellers/`. The exact file sets
(historical/ssp585 ocean chunks, OI climatology, CMIP bias, meltMIP calibration,
IMBIE2 basins, ISMIP grid, topography) are enumerated in the `OCEAN_FILES` and
`CALIBRATION_FILES` dicts at the top of `scripts/download_forcing.py` — that file
is the authoritative manifest.

No `GLOBUS_LOCAL_ENDPOINT`? The script prints the remote/local paths so you can do
the transfer by hand in the Globus web app instead.

> **Caveat (known wrinkle):** `download_forcing.py` writes to `data/forcing/...`
> using the `share_with_modellers` naming (`*_Oyr_*_ismip8km_60m_*.nc`), whereas
> the runtime expects the `ISMIP7/AIS/...` layout in §2b. Until these are unified,
> the practical path is to mirror the full `ISMIP7/AIS/` protocol tree directly
> (Globus web app) and set `ISMIP7_DATA_ROOT` to it.

### 2b. The runtime tree `ISMIP7/AIS/` (what `forcing.py` reads)

`icepack2_tools/forcing.py` resolves data as follows (override the root with
`ISMIP7_DATA_ROOT`):

```
ISMIP7/AIS/
  <ESM>/<scenario>/SDBN1-8000m/<var>/<version>/      # atmosphere (acabf, acabf-anomaly, ts, tas, pr, ...)
      <var>_AIS_<ESM>_<scenario>_SDBN1-8000m_<version>_<YEAR>.nc
  <ESM>/<scenario>/ocean/<tf|thetao|so>/<version>/   # ocean thermal forcing, salinity, temp
  <ESM>/<scenario>/fracture/                          # ice-shelf collapse / lake masks
  meltMIP/OI_Climatology_ismip8km_60m_<tf|so|thetao>_extrap.nc   # CTRL climatology
  parameterisations/ocean/imbie2/                     # IMBIE2 basin numbers (per-basin K)
  parameterisations/ocean/{bfrns,meltobs,shelfmask,floatingmasks,...}
  parameterisations/fracture/
```

Defaults baked into the readers: atmosphere `version=v2` at `SDBN1-8000m`, ocean
`version=v3`. The forcing API:

- `ISMIP7Atmosphere(esm, scenario).get_smb(year, x, y, anomaly=…)`
- `ISMIP7Ocean(esm, scenario).get_thermal_forcing(...)` / `.get_salinity(...)`
- `ISMIP7Fracture(esm, scenario).get_collapse_mask(year, x, y)`
- `make_forcing_callback(atm=, ocean=, fracture=, K=…, K_per_basin_npz=…)` — bundles
  all three into the per-step callback `run_simulation` expects.

---

## 3. Build the mesh

Adaptive isotropic mesh with Ua-style grounding-zone + calving-front refinement,
sized from BedMachine geometry and MEaSUREs strain rate (needs §1 data).

```bash
cd antarctica
python scripts/mesh_antarctica.py --lc 2500 --lc-coarse 64000 --buffer-m 20000
# or equivalently via env vars (defaults match the rest of the pipeline):
ISMIP7_LC=2500 ISMIP7_LC_COARSE=64000 ISMIP7_BUFFER_M=20000 python scripts/mesh_antarctica.py
# dev mesh used by inversion_icepack2.py / diagnostic_solve.py / run_eigendec.py:
python scripts/mesh_antarctica.py --lc 8000 --lc-coarse 80000 --buffer-m 20000
# → mesh/antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.msh
#    (e.g. antarctica_64000_2500_buffered20000.msh)
# → mesh/boundary_ids_antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.json
#    (e.g. boundary_ids_antarctica_64000_2500_buffered20000.json)
```

`--lc` / `--lc-coarse` select the fine (grounding-line/calving-front) and
coarse (interior) element sizes in meters. The GL-band element sizes,
`shelf_size`, `buffer_size`, and strain-rate floor all scale with `lc/2500`;
calving-front decay lengthscales are floored at `lc` and `1.25*lc` so they
never fall below the mesh's own fine resolution.

The mesh outline is pushed `--buffer-m` / `ISMIP7_BUFFER_M` meters into the ocean before
meshing (default `20000`; pass `0` for no buffer), which lets icepack2 handle
`h=0` at the (now-interior) calving front instead of needing a
`calving_terminus` BC. Because this changes the boundary topology, both the
`.msh` and its sidecar are named after the exact `(COARSE, FINE, BUFFER_M)`
combination used to build them — building a different resolution or buffer
never collides with or silently invalidates a previous build.

**Boundary IDs.** The gmsh physical groups come out alternating
`Calving_0, Other_1, Calving_2, …`, auto-numbered `1,2,3,…`, so **odd tag =
calving, even tag = other**. `mesh_antarctica.py` automatically writes that
split to `mesh/boundary_ids_antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.json`
(read by every solver via `ISMIP7_BNDIDS`, which itself defaults using the
same `(COARSE, FINE, BUFFER_M)`-based naming — see `scripts/mesh_naming.py`).
These small JSONs are the *only*
tracked files in `mesh/`. To regenerate one for an existing mesh without
rebuilding it, run `ISMIP7_BUFFER_M=<N> python scripts/make_boundary_ids.py`
(or pass explicit `ISMIP7_MESH`/`ISMIP7_BNDIDS` paths) — it parses the mesh's
`$PhysicalNames` block and writes `{"calving":[odd…], "other":[even…]}`.

---

## 4. Invert for basal/rheology fields (`inversion_icepack2.py`)

MAP estimate of the bed friction `θ` and rheology `φ` from the diagnostic
3-field (V × Σ × τ) system, regularized, n=1→3 continuation, via `tlm_adjoint`.

```bash
cd antarctica
ISMIP7_LC=2500 mpiexec -n 12 python scripts/inversion_icepack2.py
# → mesh/inversion_icepack2_<LC>.h5   (the MAP checkpoint every forward run loads)
```

See the script header for the full set of regularization / iteration options.

---

## 5. Calibrate ocean melt (`calibrate_melt.py`)

Solves for the Burgard quadratic-mixed-slope coefficient **K** (global and
per-IMBIE2-basin) by matching integrated observed shelf melt
(Paolo/Adusumilli ≈ 865 Gt/yr). Needs §2 forcing + §4 inversion mesh.

```bash
cd antarctica
ISMIP7_LC=2500 python scripts/calibrate_melt.py
# → results/calibrated_K_per_basin_<LC>.npz
```

The control run (§6) requires this `.npz`. Projections can either use it
(`K_per_basin_npz=`) or a scalar `ISMIP7_K_MELT`.

---

## 6. Run control & projections

All forward runs go through `scripts/simulation.py` (`setup_model` +
`run_simulation`). Drivers live in `scripts/control/` and `scripts/projections/`.

```bash
cd antarctica

# Control (CTRL2015): fixed 2000–2029 SMB climatology + OI ocean climatology,
# per-basin calibrated K, melt recomputed each step from evolving geometry.
mpiexec -n 12 python scripts/control/run.py
# → results/ctrl2015_<esm>_<lc>_{final.h5, t<year>.h5, timeseries.csv}

# Core Experiment 7: SSP5-8.5 / CESM2-WACCM, 2015–2300
mpiexec -n 12 python scripts/projections/ssp585_cesm_waccm.py
# → results/ssp585_cesm2_waccm_<lc>_{final.h5, timeseries.csv}
```

Other scenario drivers in `scripts/projections/` (ssp126/ssp370 × CESM2-WACCM /
MRI-ESM2-0, plus `ocx.py`) follow the same pattern. Historical spin-up drivers are
in `scripts/historical/`; run one first to produce
`results/hist_<esm>_<lc>_final.h5`, which the projections pick up automatically via
`ISMIP7_RESTART` (otherwise they cold-start from BedMachine).

### Environment knobs (all forward runs)

| Env var | Meaning | Default |
|---------|---------|---------|
| `ISMIP7_LC` | fine mesh resolution tag (selects mesh + inversion h5) | `2500` |
| `ISMIP7_LC_COARSE` | coarse mesh tag | `64000` |
| `ISMIP7_BUFFER_M` | outline buffer (m) used to resolve the default mesh/boundary-id filenames (see §3) | `20000` |
| `ISMIP7_MESH` | override mesh `.msh` path | `mesh/antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.msh` |
| `ISMIP7_BNDIDS` | override boundary-id JSON | `mesh/boundary_ids_antarctica_<COARSE>_<LC>_buffered<BUFFER_M>.json` |
| `ISMIP7_DATA_ROOT` | ISMIP7 forcing tree root | `<repo>/ISMIP7/AIS` |
| `ISMIP7_T_END` / `ISMIP7_DT` | end year / timestep (yr) | `2300` / `1.0` |
| `ISMIP7_OUTPUT_INTERVAL` | checkpoint every N steps | `10` |
| `ISMIP7_RESTART` | restart checkpoint | `hist_<esm>_<lc>_final.h5` if present |
| `ISMIP7_K_MELT` | scalar Burgard K (projections) | `1.15e-4` (Burgard K50) |
| `ISMIP7_K_PER_BASIN_NPZ` | per-basin K file (control) | `results/calibrated_K_per_basin_<lc>.npz` |
| `ISMIP7_ESM` | ESM for control (`CESM2-WACCM`, `MRI-ESM2-0`) | `CESM2-WACCM` |
| `ISMIP7_H_CLAMP` | thickness floor (m) | `0` |
| `ISMIP7_NO_CALVING_TERMINUS` | set to drop the calving-terminus BC | _(unset)_ |

> **dt guidance** (from a dt-convergence sweep): use `ISMIP7_DT=0.1` for
> production projections; `0.25` is acceptable if 10 steps/yr is too costly for a
> 285-yr run. `dt=1.0` over/under-melts per step and resurrects clamped cells.

---

## Outputs

Per experiment in `results/`:
- `<exp>_final.h5` — final state checkpoint (Firedrake `CheckpointFile`).
- `<exp>_t<year>.h5` — intermediate checkpoints (every `OUTPUT_INTERVAL` steps).
- `<exp>_timeseries.csv` — `year, vaf_mm_sle, mass_gt` per output step.

VAF is reported in mm of sea-level equivalent; mass in Gt.

---

## Known issues (read before trusting a long run)

- **Runaway mass gain / projection crash.** The SSP5-8.5 run to 2300 currently
  **diverges**: mass and VAF grow monotonically and accelerate (mass 24M→38M Gt,
  VAF ~56.5k→96k mm by ~2175), the diagnostic solver falls back to n=1→3
  continuation repeatedly, and the momentum solve finally fails with
  `DIVERGED_LINEAR_SOLVE` near 2175. Root cause is **not** a localized solver bug —
  it's a mass-conservation artifact: the `h_clamp_init` buffer layer + an inversion
  done *with* a thickness clamp seed the buffered ocean cells, and the composite
  rheology then feeds ice back into the melted shelves (a positive feedback). The
  fix in progress is re-inverting `θ`/`φ` with `h_clamp = 0` so the initial state
  matches true BedMachine geometry. Until then, long projections are not physically
  trustworthy. See `COMPOSITE_RHEOLOGY.md` and `CLAUDE.md` for the full writeup.
- **Forcing layout mismatch** between `download_forcing.py` (`data/forcing/`,
  `share_with_modellers` naming) and the runtime tree (`ISMIP7/AIS/`). See §2.
- **Mesh/boundary-id naming is now per-`(COARSE, FINE, BUFFER_M)`.** Mesh and
  sidecar filenames are tagged with the exact resolution and outline buffer
  used to build them (`antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.msh` /
  `boundary_ids_antarctica_<COARSE>_<FINE>_buffered<BUFFER_M>.json`, see §3),
  so a sidecar can no longer silently mismatch a mesh built with a different
  resolution or buffer. If you have older meshes/sidecars built before this
  naming convention, rename them to match or rebuild via `mesh_antarctica.py`.
- **`icepack2_tools/coupled.py`** is a WIP sketch of ice↔plume coupling and
  references a `PlumeModel` that does not yet exist in this tree — not wired into
  any run.

---

## References

- Burgard et al. 2022, *The Cryosphere* — basal-melt parameterisation assessment.
- multimelt (the reference implementation): https://github.com/ClimateClara/multimelt
- ISMIP7 ocean forcing pipeline: https://github.com/ismip/ismip7-antarctic-ocean-forcing
- Greenland companion: https://github.com/dlilien/ISMIP7_Greenland_Icepack
