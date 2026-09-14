# Forward-run readiness for the ISMIP7 submission

Swept from the ISMIP discussion board (https://github.com/orgs/ismip/discussions,
45 threads) on 13 September 2026. This file is the tracked home of that sweep.
Ordered by what blocks a submission first.

## 1. Protocol changes since 12 August

**Data freeze, 3 September (#37).** The forcing is frozen unless something is
badly wrong. The freeze copy is on Source Cooperative; Globus stays the archive
of record. Only the current version of each product is kept, and version
numbers are independent per product, so an atmosphere v2 with an ocean v3 is
fine. Use the latest of each and write the versions into the README. AIS
fracture for CESM2-WACCM ssp585 may still change, having been built from the
replaced SDBN1 v2.

**Source Cooperative mirror (#40).** Products `ismip7-ais-forcing`,
`ismip7-ais-observations`, `ismip7-ais-melt-calibration`. Re-synced with Globus
on 11 September, by hand every week or two. The forcing layout is
`data/<ESM>/<scenario>/<product>/<variable>/<file>`, with no `AIS/` level and
no version directories; the version sits in the filename. Anonymous HTTPS
listing needs a browser-like User-Agent. This is the route for NOTS, which has
no Globus endpoint; `antarctica/scripts/download_mirror.py` takes it with no
AWS tooling.

**Renamed MRI atmosphere (#37).** `SDBN1-*` became `GEMB-SDBN1-*` for
MRI-ESM2-0, a name change with the data unchanged. The core experiment uses
`SDBN1` for CESM2-WACCM and `GEMB-SDBN1` for MRI-ESM2-0; `dEBM2` is the second
downscaling for perturbed ensembles. Handled: `atmosphere_product` accepts
whichever exists.

**Control experiment (#28, #15).** ctrlclim is the 2000-2029 climatology of the
last 15 years of `historical` and the first 15 of `ssp126`, per ESM. Files are
time varying so the setup matches a projection. Fracture and excess meltwater
hold at 2015 conditions; SMB-height feedback, calving and GIA stay free.

**Historical start is free (#34).** A steady initial state may be assigned to
1960 or 1975 rather than 1850, with the 1960-1989 anomaly reference driving
changes from there. This closes the 1 K cooling question at an 1850 start.

**OCX (#32, #33, #45).** Independent of CMIP by design, spin-up method is the
group's choice, no fracture forcing for historical or OCX (use the observed
front positions in `obs/`), ice mask ends 2021. Open: the AIS OCX `dacabfdz`
regression coefficients look spatially shifted (#45); the OCX gradient fields
went to v2 on 11 September, so check the version before using SMB-height
feedback there.

**Fracture masks (#29, #30, #33).** Both ESMs' SSPs (CESM ssp585 at v2.1), none
for historical or OCX. Floating ice only. The masks light up near the grounding
line on Ross and FRIS, and the focus group suggests pairing them with a stress
criterion (Lai et al. 2020). Applied as `ISMIP7_FRACTURE=mask`.

**NaN forcing outside the downscaled mask (#39)** is intended. Zero, nearest or
a large melt are all acceptable, stated in the README. Ours fills with zero.

**CESM2-WACCM ends in 2299 (#8).** The 2300 atmosphere files were removed and
the ocean stops at 2299, while a 2015-2300 run needs 2300. Handled: the reader
persists the last year on disk exactly one year past the end and reports it
once per variable; a gap inside the series stays an error.

**Melt toolbox re-release (#25).** Rerun the calibration notebook with the new
constraint datasets. No ice-model rerun.

## 2. Local forcing, audited against the mirror (13 September)

`antarctica/scripts/audit_forcing_versions.py` compares every product the
mirror publishes for the two core ESMs against `ISMIP7/AIS` (313 GB on disk).

| what | mirror | local | verdict |
|---|---|---|---|
| CESM2-WACCM atmosphere `SDBN1-8000m` acabf, acabf-anomaly (hist, ssp126/370/585) | v2 | v2 | current |
| CESM2-WACCM ocean so, tf | v3 | v3 (ssp585 also v1 leftovers) | current, delete the leftovers |
| CESM2-WACCM fracture ssp126/370/585 | v2.1 | v2 | behind, re-download |
| MRI-ESM2-0 atmosphere | `GEMB-SDBN1-8000m` v1 | `SDBN1-8000m` v1 | data current, directory renamed |
| MRI-ESM2-0 ocean so, tf | v3 | v3 | current |
| MRI-ESM2-0 fracture | v1 | v1 | current |
| `ctrl` trees, both ESMs (atmosphere and ocean, 2015-2300) | v2, v3 | absent | download for cores 9 and 10 |
| OCX atmosphere (1979-2025) and `OCX/ocean/{main,cold,warm,vary}` (1950-2025) | v1 | absent | download for core 11; `main` is the core variant |
| ocean `extras` (ct/sa bias and climatology), `grid/ocean/ISMIP7/8km-60m` | v3 | present | fine |
| `tas`, `ts`, `pr`, `mrro`, gradients for scenarios other than ssp585; 2 km atmosphere; `dEBM2` | v2, v1 | mostly absent | needed only for SMB-height feedback or a dEBM2 member |

318 of the mirror's 358 entries are absent locally and six fracture rows are
behind. Everything the current forwards read is at the freeze version.

**Staging since that audit** (the table stays the point-in-time record): the
CESM2-WACCM fracture v2.1 files and the OCX set were downloaded on 13
September, and the ssp585 CESM2-WACCM set was staged on NOTS on 14 September.
What remains open is action 2.

## 3. Output and submission

**Status.** `ISMIP7_OUTPUT=1` makes the forward accumulate the yearly flux
means and snapshot the state each year (`icepack2_tools/ismip7_output.py`, one
Firedrake checkpoint per year, written atomically).
`antarctica/scripts/write_ismip7_output.py` regrids conservatively to the 8 km
grid through a cached supermesh overlap operator, applies the request's fill
policies and units, encodes time, and writes the 21 gridded and 10 scalar files
under `AIS/<source_id>/<ism_id>/CORE/<exp>/`. A 2-year control and a 10-year
ssp585, both at 32 km, pass every content check of
`ismip7-compliance-checker`: zero naming, numerical, spatial, attribute and
consistency errors. The experiment-length checks remain.

Conventions chosen: `acabf` is the forcing SMB with the apparent-MB correction
travelling separately as `acabf_correction`, `ligroundf` is booked into the
first floating cell, `lithk` is zero where the ice mask is zero, and
`base = orog - lithk` on the grid.

What a submission needs (#5, #16, #17, #18, #19, #20, #22, #23):

- **Variable request:** `isschecker/data/ISMIP7_variable_request.csv` in
  `ismip/ISM_SimulationChecker`; the mandatory set is ISMIP6's. Fluxes
  (`acabf`, `libmassbf`, `licalvf`, `ligroundf`) need conservative regridding
  from the DG0 cells.
- **Grids:** the description files in `ISM_SimulationChecker/gdfs` define x and
  y for every supported resolution and feed CDO directly.
- **Time encoding:** `days since 1850-01-01`. State variables are snapshots at
  1 January of the following year (2015 stamped 2016-01-01), fluxes are yearly
  means stamped 1 July, the initial state is not requested, and the filename
  year range is the years actually run. A run's `t_end` is 1 January of the
  year after its last, so projections run to 2301 and the historical to 2015,
  putting the handoff checkpoint at 2015.0.
- **Names:** `<var>_AIS_<source_id>_<ism_id>_m001_<ESM>_f001_<scenario>_C0NN_<years>.nc`
  under `AIS/<source_id>/<ism_id>/<set_id>/<set_counter>/`, `ism_id` without
  dots or underscores.
- **Checker:** bounds were relaxed in July (topg to 5500 m). `licalvf` negative
  means loss. `ligroundf` is a specific mass flux booked into the last grounded
  or first floating cell, small negatives allowed. `topg` and `lithk` must not
  be masked to the evolving ice sheet, or the scalar tool's sea-level numbers
  break.
- **Scalars:** `ismip/ismip7-scalar-processing` produces the integrated scalars
  and sea-level estimates; 2D fields and scalars share the experiment folder.
- **README:** the template is the Google document from discussion #6, recording
  forcing versions, the NaN rule, the historical start year and the spin-up.
  Drafted in `ISMIP7_README_AIS_RICE_icepack2.md`, which owns those answers;
  open items are marked `[confirm]`.
- **Upload:** email ismip6 at gmail.com with the Globus id, `AIS`, group name
  and `ism_id` to receive a folder.
- **Experiment ids:** C001 CESM2-WACCM historical, C002 MRI-ESM2-0 historical,
  C003/C004 ssp370, C005/C006 ssp126, C007/C008 ssp585 (CESM first), C009/C010
  ctrl, C011 OCX. Output goes on the standard ISMIP7 grid closest to the native
  grid, 8 km here. 3D fields are requested at a few times only and we have none.
- **Deadline:** the site still describes the 30 June round; the board refers to
  the end of September 2026 for the next one, scope to confirm.

## 4. Model-side state

**Inversions.** RC and Budd MAPs exist on the Úa-preset mesh
(`inversion_icepack2_{rc,budd}_n3_dg0_logvelnet_ua2000.h5`). Every Budd MAP
older than 13 September carries the shelf-friction defect and is unusable.

**Forward.** The RC control on the Úa mesh runs and holds (1 yr, resid 0). A
10-year ssp585 on that mesh took 10.5 minutes on 32 Sapphire Rapids ranks, 6 s
per 0.1-year step, so a 2015-2300 projection is about 5 node-hours and eleven
cores about 2.5 node-days.

**Still open:** the ctrlclim scenario default and the melt calibration rerun.
Delivered here: the collapse mask in the thickness update, the 2300 forcing
year, the `GEMB-SDBN1` path, the forcing-version audit, and the output writer.

## 5. Actions, in order

1. Run a full-length experiment through the writer, and settle the `[confirm]`
   items in the submission README draft with the group.
2. Pull the two `ctrl` trees from the mirror for cores 9 and 10. Re-run
   `audit_forcing_versions.py` before the production matrix and cite it in the
   README.
3. Set `ISMIP7_CLIM_SCENARIO=ssp126`, or read the provided `ctrl` trees.
4. Rerun the melt calibration notebook with the July toolbox. The Úa-mesh
   ssp585 hit `libmassbffl` of -0.0117 kg m-2 s-1 on small grounding-zone
   cells, past the request's -0.008 bound, which makes this concrete.
5. Optional: the stress criterion (Lai et al. 2020) alongside the collapse
   mask.
6. One forced projection end to end through the writer and the checker, then
   the matrix.
