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
| OCX atmosphere `OCX/RACMO2.3p2-ERA/SDBN1-8000m` (acabf, tas, ts, pr, mrro, gradients, 1979-2025) and `OCX/ocean/{main,cold,warm,vary}` (1950-2025) | v1 | absent | download for core 11; `main` is the core variant |
| ocean `extras` (ct/sa bias and climatology), `grid/ocean/ISMIP7/8km-60m` | v3 | present | fine |
| `tas`, `ts`, `pr`, `mrro`, gradients for scenarios other than ssp585; 2 km atmosphere; `dEBM2` | v2, v1 | mostly absent | needed only for SMB-height feedback or a dEBM2 member |

318 of the mirror's 358 entries are absent locally and six fracture rows are
behind. Everything the current forwards read is at the freeze version.

**Staging since that audit** (the table stays the point-in-time record): the
CESM2-WACCM fracture v2.1 files and the OCX set were downloaded on 13
September, and the ssp585 CESM2-WACCM set was staged on NOTS on 14 September.
On 14 September the two `ctrl` trees came down as well, 8 km atmosphere and
ocean for both ESMs, 77 GB, every transfer size-checked. A re-run of the audit
that afternoon reads 440 mirror entries, 64 of them present locally at the
mirror's own version and none behind. Every `ctrl` row for
`SDBN1-8000m`, `GEMB-SDBN1-8000m` and `ocean` (so, tf, thetao) is current, so
cores 9 and 10 have their forcing. The 376 absent entries are the 2 km
atmospheres, `dEBM2`, and the per-scenario fields that only an SMB-height
feedback or a perturbed member reads.

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

The 13 September Budd MAP is unusable as well, for a narrower reason. It was
inverted while the shelf gate still multiplied through by the grounded
indicator `He`, and the shipped gate is height above flotation alone. Measured
with `check_budd_map.py --forward` on 14 September, a diagnostic re-solve under
the shipped law reproduces that MAP's own velocity to a relative L2 distance of
0.678, where a self-consistent MAP reproduces itself to about 1e-9. The census
on the same MAP puts the old sign gate at 13 647 of 103 233 floating cells
carrying friction at the cap, the `He` form at 8 781 and the shipped gate at 0,
with 141 549 of 213 525 cells inside the `He` band. Re-inversion under the
shipped law runs as NOTS 1390416; the superseded file is kept as
`inversion_icepack2_budd_n3_dg0_logvelnet_ua2000_hegate.h5`. The RC MAP is
unaffected, since regularized Coulomb never carried the gate.

**Forward.** The RC control on the Úa mesh runs and holds (1 yr, resid 0). A
10-year CESM2-WACCM ssp585 on that mesh (NOTS job 1368723) took 10.5 minutes on
32 Sapphire Rapids ranks, 6 s per 0.1-year step, so a 2015-2300 projection is
about 5 node-hours and eleven cores about 2.5 node-days.

**Melt calibration, rerun 14 September.** The re-released calibration product
combines Paolo (2023), Davison (2023) and Adusumilli (2020), and its integrated
target is 1067.4 Gt/yr against the 865.0 Gt/yr of the Paolo plus Adusumilli
table the old calibration used. Both tables went through `calibrate_melt.py` on
the same Úa mesh, so the comparison isolates the observations:

| observations | integrated target | K* | melt at K* |
|---|---|---|---|
| Paolo, Adusumilli | 865.0 Gt/yr | 4.347e-5 | 677.0 Gt/yr |
| Paolo, Davison, Adusumilli | 1067.4 Gt/yr | 4.700e-5 | 732.0 Gt/yr |

As summary statistics of the fit, K* rises 8% and the total-match K rises with
the target by 23%, 5.553e-5 to 6.853e-5. The old table's K* on this mesh sits
within 2% of the 2500 m result that predates it, so the mesh is not what moved.
`calibrate_melt.py` takes the newer table by default, `ISMIP7_MELT_OBS_CSV`
names either, and the saved npz records which one produced it.

A forward reads the per-basin `K_basin`, each basin fitted to its own
observation, so the melt it applies moves basin by basin. The ratio of new to
old K_b spans 0.86 to 3.58 with a median near 1.15:

| basin | old K_b | new K_b |
|---|---|---|
| 1, Antarctic Peninsula fringe | 2.100e-5 | 7.515e-5 |
| 7 | 3.724e-5 | 6.314e-5 |
| 9, Amundsen | 1.464e-4 | 1.935e-4 |
| 13, the one that falls | 7.922e-5 | 6.850e-5 |

**Still open:** nothing from the August list. Delivered here: the `ssp126`
ctrlclim default (`icepack2_tools/climatology.py`), the collapse mask in the
thickness update, the 2300 forcing year, the `GEMB-SDBN1` path, the
forcing-version audit, the output writer, and the melt calibration above.

## 5. Actions, in order

1. Drive the full-length ssp585 (NOTS 1390452, 2015 to 2301 on the Úa mesh with
   `ISMIP7_OUTPUT=1`) through the writer and the compliance checker. It is the
   first run at experiment length, so it is what clears the checker's remaining
   length checks. Record which K calibration it read: a run picks up whichever
   `calibrated_K_per_basin_*.npz` is staged when it starts. The job is held in
   the queue while the slope decision of action 5 is open, and
   `scontrol release 1390452` starts it.
2. Settle the `[confirm]` items in the submission README draft with the group.
3. Carry the Budd re-inversion to a MAP that passes `check_budd_map.py
   --forward`, then run that check on every MAP the matrix will use.
4. Re-run `audit_forcing_versions.py` immediately before the production matrix
   and cite it in the README. The `ctrl` pull for cores 9 and 10 is done, and
   the mirror is re-synced with Globus by hand every week or two, so the freeze
   versions can still move under a long campaign.
5. Settle where the `libmassbffl` bound violation comes from. The request's AIS
   minimum is -0.008 kg m-2 s-1, which is 275.3 m/yr of ice, and the 10-year Úa
   ssp585 of job 1368723 reached -0.0117, or 402.6 m/yr.

   The leading explanation, now measured by `check_melt_bound.py`: the
   calibration caps the draft slope `sin(alpha)` at 5e-3 and the forward
   applies no cap, so the melt the forward applies is a different field from
   the melt the per-basin K was fitted against. The script runs both slope
   operators, capped and uncapped, on the Úa mesh at the reference geometry
   over 1 512 899 km2 of floating ice, with the per-basin K on disk:

   | slope | max, m/yr | area mean, m/yr | integrated, Gt/yr | nodes past the bound |
   |---|---|---|---|---|
   | calibration operator, capped, what K was fitted to | 71.1 | 0.77 | 1067 | 0 |
   | calibration operator, uncapped | 1804.9 | 4.18 | 5803 | 421 |
   | forward operator, uncapped, as the forward runs today | 1364.8 | 3.09 | 4293 | 298 |
   | forward operator, capped | 71.1 | 0.74 | 1028 | 0 |

   The last row is the one that settles it. Capping the forward's own operator
   integrates to 1028 Gt/yr against the 1067.4 Gt/yr the K was fitted to
   reproduce, within 4%, and leaves nothing past the bound. Running uncapped
   integrates four times the target. The two operators differ because the
   calibration projects `grad(draft)` from a CG1 geometry while the forward
   differentiates a `cg1_lift` of a DG0 draft, which is smoother; the cap
   removes that difference, since both then sit at 5e-3 almost everywhere.

   The 10-year run's own budget, 1860 Gt/yr, sits below the uncapped reference
   value, so its melt-receiving mask and its evolved geometry account for part
   of the difference as well. A node past the bound has a median area of
   6.33 km2 against 64 km2 for an 8 km pixel.

   The slope gap leaves two candidates open, and the full-length run still has
   to settle them: the evolved geometry with its warmer projected thermal
   forcing, and the bookkeeping, since `book_advance` books the melt REQUESTED
   of a step while a nearly ice-free floating cell can only lose what it holds,
   after which the request's `no_floating_ice` fill policy reports that cell's
   rate for the whole 8 km pixel.

   Closing the gap between calibration and forward is the next step. Which side
   moves is a science decision, since the cap is tied to the unsettled upstream
   local-slope question. Until it is made, `load_K_per_basin` warns once per
   run when the K file it reads records the cap it was fitted against.
6. Optional: read the provided `ctrl` trees in place of the `ssp126`
   reference-climate pool.
7. Optional: the stress criterion (Lai et al. 2020) alongside the collapse
   mask.
