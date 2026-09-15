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

   `check_melt_bound.py` measures the slope side of it. The calibration caps the
   draft slope `sin(alpha)` at 5e-3 and the forward applies no cap, so the melt
   the forward applies is a different field from the melt the per-basin K was
   fitted against. The script evaluates two halves at the reference geometry
   with `calibrated_K_per_basin_2000.npz` on the Úa 2 km mesh, each capped and
   uncapped. The calibration half reproduces `calibrate_melt.py` on CG1 nodes,
   with BedMachine's raster surface and mask and the cap on the nodal slope.
   The forward half reproduces the forward on DG0 cells, with the surface from
   flotation, `forcing.compute_sin_alpha`'s cell slope, forcing at each cell's
   own draft, the callback's `haf <= 0` floating test and the cap on the cell
   slope. The calibration half floats 1 512 899 km2 over 47 288 nodes and the
   forward half 1 631 466 km2 over 85 820 cells, and they differ in mask and
   quadrature as well, so their totals compare in magnitude.

   | slope | max, m/yr | p99, m/yr | area mean, m/yr | integrated, Gt/yr | past the bound |
   |---|---|---|---|---|---|
   | calibration half, capped, what K was fitted to | 71.1 | 22.2 | 0.77 | 1067 | 0 |
   | calibration half, uncapped | 1804.9 | 256.3 | 4.18 | 5803 | 421 nodes, 2920.6 km2 |
   | forward half, uncapped, as the forward runs today | 1144.1 | 59.0 | 1.16 | 1732 | 97 cells, 388.2 km2 |
   | forward half, capped | 57.8 | 13.0 | 0.43 | 646 | 0 |

   The same script with `calibrated_K_per_basin_2500.npz`, the coefficient file
   the 10-year run of job 1368723 read, everything else unchanged, gives the
   like-for-like comparison against that run:

   | slope, 2500 file | max, m/yr | area mean, m/yr | integrated, Gt/yr | past the bound |
   |---|---|---|---|---|
   | calibration half, capped, what K was fitted to | 54.0 | 0.61 | 841 | 0 |
   | calibration half, uncapped | 1522.7 | 3.34 | 4634 | 320 nodes, 2193.7 km2 |
   | forward half, uncapped, as the forward runs today | 869.8 | 0.92 | 1380 | 63 cells, 248.6 km2 |
   | forward half, capped | 43.9 | 0.34 | 510 | 0 |

   The calibration half reproduces the target its own K was fitted to: 1067
   against 1067.4 Gt/yr for the 2000 file, and 841 against 865 Gt/yr for the
   2500 file, within 3%. That is the internal consistency check. With the same
   K, the forward applies 1.62 times the target with the 2000 file, 1732
   against 1067, and 1.64 times with the 2500 file, 1380 against 841. The ratio
   is the durable result, stable across both coefficient files.

   The 10-year run booked 1860 Gt/yr with the 2500 file, above the 1380 Gt/yr
   its forward half applies at the reference state, since that run carries
   warmer ssp585 thermal forcing over evolving geometry where the script holds
   the OI climatology at the initial state. That gap is a consistent residual.
   An earlier reading of this section put the uncapped forward half within 7%
   of the run; it compared the 2000 file's 1732 Gt/yr with a run on the 2500
   file and does not hold.

   Capping the forward's own slope gives 646 Gt/yr with the 2000 file, 39%
   under the target. Neither convention on its own reconciles the halves,
   which also differ in floating area, mask and quadrature.

   An earlier form of the script lifted the forward's slope onto CG1 nodes and
   melted it with CG1 forcing and the raster mask. Its forward rows, 4293 Gt/yr
   and 298 nodes past the bound uncapped and 1028 Gt/yr capped, reproduced
   neither half and are superseded, and with them the reading that capping the
   forward lands within 4% of the target.

   Both mechanisms for the bound violation are visible at the reference state.
   Uncapped, 97 forward cells melt past the bound over 388.2 km2, so the
   parameterisation itself exceeds it there. Their median area is 3.61 km2
   against 64 km2 for an 8 km pixel, so the gridded value comes from a small
   hot cell filling its pixel under the request's `no_floating_ice` fill
   policy. The full-length run still has to settle two further contributions:
   the evolved geometry with its warmer projected thermal forcing, and the
   bookkeeping, since `book_advance` books the melt REQUESTED of a step while a
   nearly ice-free floating cell can only lose what it holds.

   Closing the gap between calibration and forward is the next step. The clean
   route is to recalibrate K through the forward's own melt path, cell by cell
   with its own floating mask, under whichever slope convention is chosen.
   Choosing the convention is a science decision, since the cap is tied to the
   unsettled upstream local-slope question. Until it is made, `load_K_per_basin`
   warns once per run when the K file it reads records the cap it was fitted
   against.
6. Optional: read the provided `ctrl` trees in place of the `ssp126`
   reference-climate pool.
7. Optional: the stress criterion (Lai et al. 2020) alongside the collapse
   mask.
