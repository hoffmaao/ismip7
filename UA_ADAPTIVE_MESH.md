# Úa-style adaptive remeshing (branch `ua-adapt-mesh`)

Port of the mesh adaptation Úa runs between run-steps, read from UaSource
(`UaMain/AdaptMesh.m`, `NewDesiredEleSizesAndElementsToRefineOrCoarsen2.m`,
`Error2EleSize.m`, `GlobalRemeshing.m`, `MapFbetweenMeshes.m`, Sep 2026).

## What was ported, and what was not

**Not ported: `GLmorphing`.** Úa's repository carries a mesh-deformation
scheme (`GLmorphing.m` / `GLmorphingInit.m`): every node is frozen at fixed
barycentric coordinates inside a constrained Delaunay triangulation of the
domain boundary plus the grounding-line polyline, and nodes follow the
grounding line as those vertices move. It has no caller in `Ua.m` or `Ua2D.m`,
its hook in `GlobalRemeshing.m` is commented out with "broken anyhow", and the
362-page Compendium never mentions it. What Úa runs is `AdaptMesh`, and that
is what this branch implements. The morphing idea is ~100 lines on top of
`icepack2_tools.adapt_mesh.grounding_line_points` if it is ever wanted.

## The scheme (Úa names in brackets)

1. **Desired element size** at the nodes of the current mesh
   [`EleSizeDesired`]. Start at `MeshSizeMax`. For each enabled criterion
   [`ExplicitMeshRefinementCriteria`] compute a nodal error proxy `e` and map
   it [`Error2EleSize`]: `h = hMin + (e0/(e+e0))^(1/p) (hMax - hMin)`, with
   `e0` the criterion's `Scale`. Minimum over criteria; if none fired, use
   `MeshSize`. Relax toward the current size with `W = 0.5` and clip the ratio
   of change to `[1/5, 5]`. THEN the absolute bands: nodes within `d_i` of the
   grounding line get `min(h, s_i)` per row of `MeshAdapt.GLrange`, floored at
   `MeshSizeMin`; likewise `CFrange` at the calving front.
2. **Global remeshing** [`explicit:global`, gmsh]: the nodal size field on the
   old mesh becomes gmsh's background scalar view, the domain outline and
   physical groups are rebuilt exactly as `mesh_antarctica.py` builds them,
   and the mesh is regenerated. Then Úa's element-count control: up to four
   rescalings of `MeshSizeMin` so the count lands within
   `[LowerLimitFactor, UpperLimitFactor] * MaxNumberOfElements`.
3. **Transfer** [`MapFbetweenMeshes`]: point-evaluation interpolation old to
   new (Úa: FE shape functions), thickness `ThickMin` outside the old mesh,
   the bed re-sampled from BedMachine on the new mesh (Úa's `DefineGeometry`
   route, with the sampling method the MAP recorded), surface by flotation.
   Frozen anchors (`a_ref_mb`, `N_ref`, `C_w0`, `H_init`, level set) move too.

Criteria available, exactly Úa's list: `effective strain rates`,
`effective strain rates gradient`, `flotation` (Úa's `DiracDelta`
`0.5 k sech^2(k (h - h_f))`), `thickness gradient`, `upper surface gradient`,
`lower surface gradient`, `|dhdt|`, `dhdt gradient`.

Two deliberate choices where Úa is nodal and this model is DG0:

- `EleSizeCurrent` is Úa's `sqrt(mean adjacent element area)`, kept even
  though it is 0.66 of an edge length, because Úa's relaxation and ratio
  limits are calibrated against it.
- Thickness transfer: Úa moves its NODAL surface; for a DG0 thickness the
  analogue that keeps mass is the conservative supermesh projection
  (`ISMIP7_ADAPT_TRANSFER=project`, the recommended setting; it must run on
  one rank, see below). The surface route lost 3.8% of the volume on a 130 km
  interior; the projection kept it to 0.01%. Either way the transfer prints
  the volume change and the mean front thickness before and after.

## Configuration (`ISMIP7_ADAPT_*`, defaults = Úa's `Ua2D_DefaultParameters`)

| variable | Úa `CtrlVar` | default |
|---|---|---|
| `ISMIP7_ADAPT_MESH_SIZE` | `MeshSize` | 10 km |
| `ISMIP7_ADAPT_MESH_SIZE_MIN` / `_MAX` | `MeshSizeMin` / `MeshSizeMax` | 1 km / 10 km |
| `ISMIP7_ADAPT_GL_RANGE="5000:2000,1000:500"` | `MeshAdapt.GLrange=[5000 2000; 1000 500]` | none |
| `ISMIP7_ADAPT_CF_RANGE` | `MeshAdapt.CFrange` | none |
| `ISMIP7_ADAPT_CRITERIA="effective strain rates:0.01,..."` | `ExplicitMeshRefinementCriteria(I).Name/Scale[/p/EleMin/EleMax]` | none |
| `ISMIP7_ADAPT_RELAXATION_W` | `W` | 0.5 |
| `ISMIP7_ADAPT_MAX_RATIO_CHANGE` / `_MIN_` | `Max/MinRatioOfChangeInEleSizeDuringAdaptMeshing` | 5 / 0.2 |
| `ISMIP7_ADAPT_MAX_ELEMENTS` | `MaxNumberOfElements` (0 = off) | 0 |
| `ISMIP7_ADAPT_THICK_MIN` | `ThickMin` | 1 m |
| `ISMIP7_ADAPT_DIRAC_WIDTH` | `RefineDiracDeltaWidth` | 100 m |
| `ISMIP7_ADAPT_TRANSFER` | (none) | `interpolate` |
| `ISMIP7_ADAPT_GEOMETRY` | `MapOldToNew.Transient.Geometry` | `bh-FROM-sBS` |
| `ISMIP7_ADAPT_KEEP_CURRENT=1` | (none; null test at the current sizes) | off |

### The Úa preset: `ISMIP7_ADAPT_PRESET=ua`

Sets the sizes Úa itself uses for Antarctica, so the adaptation is consistent
with a Úa run rather than with this repo's initial meshes. Provenance per value:

| setting | value | source |
|---|---|---|
| `MeshSizeMax` (interior) | 180 km | Úa-FESOM pan-Antarctic mesh, GMD 18 (2025): "up to 180 km in the interior" |
| `MeshSizeMin` (grounding line) | 2 km | same paper: "adaptive refinement down to 2 km at the grounding line" |
| `MeshSize` (fallback) | 90 km | PIG-TWG example: `MeshSize = MeshSizeMax/2` |
| ice-shelf size | 4 km | PIG-TWG: `MeshSizeIceShelves = MeshSizeMax/5`; pan-Antarctic "4 km" |
| low ground (`s < 1500 m`) | 36 km | PIG-TWG: `EleSizeIndicator(s<1500) = MeshSizeMax/5`, scaled to the 180 km max |
| `effective strain rates` criterion | Scale 0.001, floor 4 km | PIG-TWG Scale; pan-Antarctic "4 km in regions of high strain rate" |
| `GLrange` | 10 km: 4 km, 5 km: 2 km | INFERRED pan-Antarctic form of MISMIP+'s `[20000 5000; 10000 2000; 5000 500]` and the "2 km at the GL" statement |
| `MaxNumberOfElements` | 250,000 | pan-Antarctic Úa: "250 000 elements" (linear, `TriNodes=3`, as here) |
| initial adaptation | iterated up to 5x | PIG-TWG: `AdaptMeshMaxIterations=5` |
| remesh interval | every step in Úa | here `--adapt-every` in years; 1 yr is the practical floor |

Override any single value with its `ISMIP7_ADAPT_*` variable; the preset only
fills what is unset. Úa's elements in these setups are linear (`TriNodes=3`),
so a 4 km Úa element and a 4 km cell here resolve alike.

Timing knobs live in the orchestrator (`run_adaptive.py`): `--adapt-every`
[`AdaptMeshTimeInterval`], `--initial-iterations`
[`AdaptMeshInitial` + `AdaptMeshMaxIterations`].

## Pieces

- `icepack_tools/adapt_mesh.py` (the shared package, next to the level set) -
  the scheme itself: criteria, `Error2EleSize`, relaxation and ratio limits,
  bands, the PIG-TWG rules, `remesh_global` with a project-supplied
  `build_geometry()` and Úa's element-count control, and the transfer helpers
  (`cross_mesh_transfer`, `preserve_front`, `physical_divergence`,
  `surface_route_thickness`, `rebuild_reference_pressure`). Nothing in it
  knows about Antarctica. Tests: `icepack_tools/test/adapt_mesh_test.py`.
  Import it before assembling any form: it reaches `icepack2` through
  `icepack_tools.constants`, and Irksome refuses to load afterwards.
- `icepack2_tools/adapt_mesh.py` - what ISMIP7 adds: the Antarctic domain
  builder for the remesh and `transfer_state`, which knows this model's
  checkpoint (its fields, frozen references and restart attributes).
- `antarctica/scripts/adapt_mesh.py CHK --out-checkpoint NEW [--rebuild-aref]`
  - one adaptation of a forward checkpoint; writes the new `.msh`, its
  boundary-id sidecar and a restartable checkpoint.
- `antarctica/scripts/run_adaptive.py` - the segment loop: run the ordinary
  driver to the next adaptation time, adapt, restart.
- `simulation.py`: a restart whose checkpoint carries `adapted_initial=1` and
  no `a_ref_mb` rebuilds the apparent-mass-balance reference on the new mesh.

## Geometry route

Úa's default `MapOldToNew.Transient.Geometry = "bh-FROM-sBS"` moves the
**surface** and derives the thickness from it and the re-sampled bed
(`h = min(s - b, s rho_w / (rho_w - rho_i))`). The alternative
`"bs-FROM-hBS"` moves the thickness. Measured on the 32 km control refined to
8 km bands, the thickness route put 32 km cell-mean thicknesses onto a finely
re-sampled bed and floated them over every trough (outflux 782 -> 14,770 Gt/yr
at the first step), so the default is Úa's. The DG0 analogue of Úa's nodal
surface is the volume-preserving CG1 lift of the cell surface, point-evaluated
at the new centroids. On the FIRST run-step of a transient run Úa takes all
geometry from data; `--rebuild-aref` does the same (bed and thickness from
BedMachine, `H_init = H`).

## Two things a nodal model never had to face

**The frozen mass-balance reference does not transfer.** `a_ref_mb` cancels
the *old* mesh's discrete flux divergence spike by spike (about 1000 m/yr at
the Pine Island grounding zone). Interpolated onto another mesh those spikes
become misplaced sources: a same-resolution remesh blew up within a year
(outflux 4,371 -> 53,185 Gt/yr). So the transfer carries the *physical*
divergence instead, `P = flux/area - a_ref` (conservatively projected, saved
as `phys_div`), and the forward rebuilds `a_ref = flux_new/area - P` with its
own upwind operator on the new mesh, exactly as it built the t=0 reference.

**Budd's effective-pressure reference is frozen to the original geometry.**
`N_ref` is the t=0 effective pressure and the law reads `N_hat = N_eff/N_ref`
(cap 3, zero afloat). Re-sampling the bed on a new mesh changes `N_eff` cell
by cell while an interpolated `N_ref` does not follow. Measured on a
same-resolution null remesh, the grounded-cell ratio spread went from
0.91 to 1.48 (5th to 95th percentile) to 0.36 to 2.97, with 10% of cells
losing their reference altogether, which the law reads as triple friction;
on the 8 km-band case 26% did. The transfer therefore carries the RATIO the
run had and rebuilds `N_ref = N_eff_new / ratio` on the new geometry, in the
law's own units. Regularized Coulomb has no such reference (its cap uses the
live `N`), so RC forwards are immune to this one.

**The DG0 front thickness smears under any interpolation.** The last cell at
the calving front is thin; a new front cell whose centroid lands in an old
interior cell inherits a thick value, and the terminus traction goes as h^2
(the geometry.py docstring measured a 4.7x outflux multiplier from the same
effect in the CG1 lift). The transfer prints the mean front thickness before
and after; `ISMIP7_ADAPT_TRANSFER=project` (supermesh) bounds it by overlap
weighting; `ISMIP7_ADAPT_KEEP_CURRENT=1` remeshes at the current sizes to
measure the transfer's own cost with no refinement.

## Found on the way: Budd shelf friction was a sign test on roundoff

The identity transfer (`--no-remesh`) reproduced every checkpoint field to
machine precision, yet two identity checkpoints whose `surface` differed by
2e-13 m gave 672 and 1459 Gt/yr for the same year. Neither number was the
transfer's fault. In `build_rc_residual(fric_law="budd")` the shelf test was
`conditional(gt(N, 0), N_hat, 0)`: on a floating cell the surface *is* the
flotation branch, so `N = max(p_I - p_W, 0)` is a roundoff residue of either
sign, and a positive residue passed the test with `N_ref` equally tiny, where
the delta floor `nhat_floor * p_I / N_ref` lifts it to the cap. Result: triple
Weertman friction on whichever shelf cells happened to round positive (445 of
3515 in the 32 km control; the 2e-13 change flipped all of them). The shared
`icepack_tools.friction.basal_stress` already had the fix, an `He` gate
(height above flotation, 0 well below flotation whatever the roundoff does);
it is now ported here. Verified: the two checkpoints give outflux 1478 vs 1475
Gt/yr and VAF 57638.204 vs 57638.212 mm SLE. The outflux numbers quoted above
(782, and the identity test's 672) were measured with the old test in place,
as was every Budd forward and Budd MAP to date; regularized Coulomb has a
continuous `tau_cap` and never had the problem.

**Correction (same day):** the He gate alone was not enough. `He` is a smooth
function of height above flotation, so a cell floating by a few metres sits
inside the He band and still received `He * nhat_cap` from a roundoff-positive
`N` (133 of the 3791 floating cells of the 32 km MAP). The gate is now HAF > 0
itself (`dual_friction.budd_nhat`), which for grounded ice is the same
statement as N > 0 (`N = rho_I g HAF` when `s = b + H`) and on the shelf is a
real negative number instead of a cancelling difference. Census script:
`antarctica/scripts/check_budd_map.py MAP [--forward]` (old gate 418 cells at
the cap, He-only 133, production 0; and the old MAP's saved velocity is not
reproduced by the fixed law, rel L2 = 0.91, so Budd MAPs must be re-inverted).
