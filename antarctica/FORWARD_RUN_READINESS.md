# Forward-run readiness for the ISMIP7 submission

Checked against the ISMIP discussion board (https://github.com/orgs/ismip/discussions,
45 threads, all read) on 13 September 2026. The previous sweep is the "Upstream
ISMIP7 protocol status (checked Aug 12)" section of the local `CLAUDE.md`
working notes, which this repository does not track; this file supersedes it
for anything about running and submitting the forward matrix, and is the
tracked home of the sweep from here on.
Items are ordered by what blocks a submission first.

## 1. What the protocol now says (new or changed since 12 August)

- **Data freeze, 3 September (#37).** The forcing is frozen "unless someone
  finds something really wrong". The freeze copy is the one on Source
  Cooperative; Globus stays the archive of record. Only the current version
  of each product is kept. Version numbers are independent per product, so
  mixing an atmosphere v2 with an ocean v3 is fine; use the latest of each at
  submission time and **write the versions you used into the README**. One
  caveat: AIS fracture for CESM2-WACCM ssp585 may still change (it was built
  from SDBN1 v2, which was replaced over CESM topography problems).
- **Source Cooperative mirror (#40).** `aws s3 ls --no-sign-request
  --endpoint-url https://data.source.coop s3://ismip/`; products
  `ismip7-ais-forcing`, `ismip7-ais-observations`, `ismip7-ais-melt-calibration`.
  Re-synced with Globus on 11 September; the maintainer syncs by hand every
  week or two. The layout inside the forcing product is
  `data/<ESM>/<scenario>/<product>/<variable>/<file>` (no `AIS/` level and
  no version directories; the version is in the filename, see section 2).
  Anonymous HTTPS listing works with a browser-like User-Agent (Python's
  default one is refused). This is the practical route for NOTS, which has no
  Globus endpoint: install `awscli` in the venv and `aws s3 sync
  --no-sign-request` the scenario trees straight to `/projects`.
- **Renamed MRI atmosphere directories (#37).** `SDBN1-*` became
  `GEMB-SDBN1-*` for MRI-ESM2-0 (a name change only, data unchanged: v2 users
  need not rerun). The core experiment uses `SDBN1` for CESM2-WACCM and
  `GEMB-SDBN1` for MRI-ESM2-0; `dEBM2` is the second downscaling for
  perturbed ensembles. Our reader looks for `SDBN1-8000m`; it must accept the
  new name for MRI.
- **Control experiment (#28, #15).** ctrlclim is the 2000-2029 climatology of
  the last 15 years of `historical` and the first 15 of `ssp126`, per ESM;
  files are time-varying (repeated entries) so the setup is identical to a
  projection; `acabf-anomaly`/`tas-anomaly` for `ctrl` were being added
  (ctrl `mrro` landed 11 September); fracture and excess meltwater held at
  2015 conditions; SMB-height feedback, calving and GIA free. Our
  `ISMIP7_CLIM_SCENARIO` still defaults to ssp585: change it to ssp126 or
  read the provided `ctrl` tree (absent locally, see section 2).
- **Historical start is free (#34).** A steady initial state may be assigned
  to 1960 or 1975 rather than 1850; the 1960-1989 anomaly reference then
  drives changes relative to that state. This closes the question we had
  flagged about the ~1 K cooling at an 1850 start.
- **OCX (#32, #33, #45).** Independent of CMIP by design; spin-up method is
  the group's choice; **no fracture forcing for historical or OCX** (use the
  observed front positions in `obs/`, or calibrate); the OCX ice mask ends in
  2021. Open and unanswered: the AIS OCX `dacabfdz` regression coefficients
  look spatially shifted (#45, 10 September; the OCX gradient fields went to
  v2 on 11 September, so check the version before using SMB-height feedback
  in OCX).
- **Fracture masks (#29, #30, #33).** Available for both ESMs' SSPs (CESM
  ssp585 at v2.1, ours is v2), none for historical/OCX. Apply to floating ice
  only; the masks light up near the grounding line on Ross and FRIS, and the
  focus group suggests combining them with a stress criterion (Lai et al.
  2020). Still not applied in our thickness update.
- **NaN forcing outside the downscaled mask (#39)** is intended; filling with
  zero, nearest or a large melt are all acceptable, to be stated in the
  README. Ours fills with zero (`forcing.py`, `nan_to_num(nan=0.0)`).
- **CESM2-WACCM ends in 2299 (#8).** The year-2300 atmosphere files were
  removed and the ocean forcing stops at 2299. A 2015-2300 run needs the
  2300 forcing year: done, the reader persists the last year on disk exactly
  one year past the end of the series and reports it once per variable, while
  a gap inside the series stays an error. (The 2290-2299 mean was the other
  option; the steering committee had not decided.)
- **Melt toolbox re-release (#25)** is unchanged since August: rerun the
  calibration notebook with the new constraint datasets; no ice-model rerun.

## 2. Forcing we hold locally, audited against the mirror (13 September)

`antarctica/scripts/audit_forcing_versions.py` lists every product the mirror
publishes for the two core ESMs and compares versions with `ISMIP7/AIS` (313
GB on disk). The mirror keeps no version directories; the version is in the
filename (`_v2_`, fracture `-v2.1.nc`).

| what | mirror (current) | local | verdict |
|---|---|---|---|
| CESM2-WACCM atmosphere `SDBN1-8000m` acabf, acabf-anomaly (hist, ssp126/370/585) | v2 | v2 | current |
| CESM2-WACCM ocean so, tf (hist, ssp126/370/585) | v3 | v3 (ssp585 also v1 leftovers, thetao v3) | current; delete the v1 leftovers |
| CESM2-WACCM fracture ssp126/370/585 | **v2.1** | v2 | **behind: re-download** (collapse mask, lake properties, excess melt) |
| MRI-ESM2-0 atmosphere (hist, ssp126/370/585) | `GEMB-SDBN1-8000m` v1 | `SDBN1-8000m` v1 | data current, directory renamed: the reader must accept `GEMB-SDBN1` |
| MRI-ESM2-0 ocean so, tf | v3 | v3 | current |
| MRI-ESM2-0 fracture | v1 | v1 | current |
| `ctrl` trees (both ESMs): atmosphere acabf, acabf-anomaly, ts, ts-anomaly, mrro, mrro-anomaly, gradients, 2015-2300; ocean so, tf, thetao | v2 / v3 | absent | **download for cores 9 and 10** (the protocol's own ctrlclim, instead of building it from historical + ssp126) |
| OCX: `OCX/RACMO2.3p2-ERA/SDBN1-8000m` (acabf, tas, ts, pr, mrro, gradients, 1979-2025) and `OCX/ocean/{main,cold,warm,vary}` (so, tf, thetao 1950-2025, v1) | v1 | absent | **download for core 11**; the ocean has four expert-judgment variants, `main` is the core one |
| ocean `extras` (ct/sa bias and climatology, v3, 1995-2024), `grid/ocean/ISMIP7/8km-60m` v3 | v3 | present (obs, grid) | fine |
| `tas`, `ts`, `pr`, `mrro`, `dacabfdz`, `dtsdz`, `dmrrodz` for scenarios other than ssp585; the 2 km atmosphere; `dEBM2` | v2 / v1 | mostly absent | not needed unless SMB-height feedback or the dEBM2 ensemble member is run |

318 of the mirror's 358 entries are absent locally; only the six fracture
rows are behind. Everything the current forwards read (SMB, SMB anomaly,
so, tf) is at the freeze version.

## 3. Output and submission

**Status (13 September, this branch):** the writer exists and passes the
compliance checker's content checks. `ISMIP7_OUTPUT=1` makes the forward
accumulate the yearly flux means and snapshot the state each year
(`icepack2_tools/ismip7_output.py`, one Firedrake checkpoint per year,
written atomically so an interrupted run cannot damage the years already
banked);
`antarctica/scripts/write_ismip7_output.py` regrids conservatively to the
8 km grid through a cached supermesh overlap operator, applies the request's
fill policies and units, encodes time, and writes the 21 gridded and 10
scalar files with the protocol names under `AIS/<source_id>/<ism_id>/CORE/<exp>/`.
On a 10-year ssp585 at 32 km the checker (`ismip7-compliance-checker`, a
Python 3.13 venv wrapped in `~/.local/bin`) reports 0 naming, numerical,
spatial, attribute and consistency errors; what remains is the experiment
length. Conventions chosen: `acabf` is the forcing SMB and the apparent-MB
correction travels separately as `acabf_correction` (not a request
variable), `ligroundf` is booked into the first floating cell, `lithk` is
zero where the ice mask is zero, and `base = orog - lithk` on the grid.
Still to do: a full-length run through it, the scalar tool cross-check, the
group's decisions on the `[confirm]` items in the submission README draft
(`ISMIP7_README_AIS_RICE_icepack2.md`), and the submission email.

What a submission needs (#5, #16, #17, #18, #19, #20, #22, #23):

- **Variable request:** `isschecker/data/ISMIP7_variable_request.csv` in
  `ismip/ISM_SimulationChecker` (the old `conventions/` path is gone); the
  mandatory set is ISMIP6's. `icepack2_tools/regrid.py` already carries the
  variable table and a point-sampling regridder onto the 8 km grid; fluxes
  (acabf, libmassbf, licalvf, ligroundf) must be regridded conservatively
  from the DG0 cells, not point-sampled.
- **Grids:** the grid description files in `ISM_SimulationChecker/gdfs`
  define x/y for every supported resolution and feed CDO directly.
- **Time encoding:** `days since 1850-01-01`; state variables are snapshots
  at 1 January of the following year (2015 output stamped 2016-01-01), flux
  variables are yearly means stamped 1 July; the initial state is not
  requested; the filename year range is the years actually run
  (`..._C007_2015-2300.nc`, historical `..._C001_<start>-2014.nc`). A run's
  `t_end` is 1 January of the year AFTER the last one it covers, since `t=Y.0`
  is 1 January of year Y: the projections run to 2301 and the historical to
  2015, which is also what puts the handoff checkpoint at 2015.0 so the
  projection's first forcing year is 2015.
- **Names:** `<var>_AIS_<source_id>_<ism_id>_m001_<ESM>_f001_<scenario>_C0NN_<years>.nc`
  under `AIS/<source_id>/<ism_id>/<set_id>/<set_counter>/`; `ism_id` without
  dots or underscores.
- **Checker:** `ISM_SimulationChecker` must pass (bounds relaxed in July:
  topg to 5500 m; `licalvf` negative means loss; `ligroundf` is a specific
  mass flux booked into the last grounded or first floating cell, small
  negatives allowed; `topg` and `lithk` must not be masked to the evolving
  ice sheet or the scalar tool's sea-level numbers break).
- **Scalars:** `ismip/ismip7-scalar-processing` produces the integrated
  scalars and sea-level estimates from the 2D fields; 2D fields and scalars
  go in the same experiment folder.
- **README:** the template is the Google document linked from
  discussion #6; it must record forcing versions, the NaN rule, the
  historical start year and the spin-up. Drafted against it in
  `ISMIP7_README_AIS_RICE_icepack2.md`, which is what gets submitted and is
  the owner of those answers; everything still open there is marked
  `[confirm]`.
- **Upload:** email ismip6 at gmail.com with the Globus id, `AIS`, the
  group name and `ism_id` to receive an upload folder.
- **Experiment ids** (ismip.org/research/ismip7, read 13 September): C001
  CESM2-WACCM historical, C002 MRI-ESM2-0 historical, C003/C004 ssp370,
  C005/C006 ssp126, C007/C008 ssp585 (CESM first, MRI second), C009/C010
  ctrl, C011 OCX. Output goes on "the standard ISMIP7 grid that is closest
  to a model native grid" (8 km for us); 3D fields are only requested at a
  few times and we have none.
- **Deadline:** the page still describes the 30 June round (C007 and C001
  through the checker); the board refers to "the deadline at the end of the
  month" (September 2026) for the next round, scope to confirm on the status
  spreadsheet.

## 4. Model-side state (this branch)

- Inversions: RC and Budd MAPs exist on the Úa-preset mesh
  (`inversion_icepack2_{rc,budd}_n3_dg0_logvelnet_ua2000.h5`); the Budd one is
  being re-inverted under the shipped friction law (job 1339328). Every Budd
  MAP older than 13 September carries the shelf-friction defect and is not
  to be used.
- Forward: the RC control on the Úa mesh runs and holds (1 yr, resid 0);
  no forced projection has been run on that mesh yet. Cost estimate from the
  inversion's solves (140 s per diagnostic solve on 32 Sapphire Rapids
  ranks, 10 steps per year): about 25 minutes per simulated year, so a
  2015-2300 projection is roughly 5 node-days; the `commons` 1-day limit
  means 5 chained links per experiment, the `long` 3-day limit 2. Eleven
  cores at that cost are about 55 node-days.
- Protocol wiring still open: ctrlclim scenario (ssp126), and the melt
  calibration rerun. Delivered on this branch: the collapse mask in the
  thickness update (`ISMIP7_FRACTURE=mask`), the 2300 forcing year, the
  `GEMB-SDBN1` path for MRI, the forcing-version audit, and the ISMIP7
  output writer (section 3).

## 5. Actions, in order

1. Output writer: done for the content checks (section 3); the README is
   drafted (`ISMIP7_README_AIS_RICE_icepack2.md`). Next, run a full-length
   experiment through it and settle the draft's `[confirm]` items with the
   group.
2. Still to pull from the mirror with
   `antarctica/scripts/download_mirror.py`: the two `ctrl` trees and the
   remainder of OCX. (CESM2-WACCM fracture v2.1 and the OCX set are staged
   locally, and 39 GB of ssp585 forcing is staged on NOTS; the downloader
   replaces the `awscli` route, needs no AWS tooling, and resumes.) Re-run
   `audit_forcing_versions.py` before the production matrix and cite its
   output in the README.
3. End-of-series rule for 2300 in `forcing.py` and the `GEMB-SDBN1` path for
   MRI: both done (`_load_year` bridges exactly one year past the end;
   `atmosphere_product` accepts either product name).
4. `ISMIP7_CLIM_SCENARIO=ssp126` default (or the provided `ctrl` trees).
5. Fracture masks on floating ice in the thickness update: done as
   `ISMIP7_FRACTURE=mask` (floating cells the mask flags are emptied and
   booked as calving; CESM2-WACCM masks re-downloaded at v2.1); the stress
   criterion (Lai et al. 2020) is not implemented and stays optional.
6. Rerun the melt calibration notebook with the July toolbox.
7. One forced projection on the Úa mesh end to end through the writer and
   the checker, then the matrix.
