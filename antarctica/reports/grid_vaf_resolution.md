# The 8 km grid's volume above flotation, by mesh resolution (issue #99)

IU Quartz, 26 September 2026. ismip7-scalars takes volume above flotation
(VAF) from the 8 km pixel means of `lithk` and `topg`, and in the 32 km p2
ssp585 run that puts its slvaf 39.4 mm above the model's own by 2300, about
10 % of that run's signal (`scalar_comparison_32km.md`, cause 2). This asks
whether a finer mesh brings the difference down to the 1 to 2 % the tool's
author reported (forum thread 19). No run finer than 32 km on Quartz has
annual output past its first year, so the check takes the 2015 state of every
mesh from 32 km to 0.5 km, changes it the same way on each, and measures the
tool's error directly.

It does not come down. Removing every ice shelf moves the tool's slvaf by
+17.6 mm at 32 km and +20.3 mm on the production layout. The 32 km run's
change by 2300, imposed on each mesh, moves it by +43.0 mm at 32 km (the run
itself reached +39.4 mm) and +55.8 mm on the production layout.

## Method

`scripts/grid_vaf_attribution.py` builds the writer's overlap operator for
a mesh, forms the pixel means the writer submits (`lithk` over the whole
pixel, `topg` over the covered part), and sets the tool's integrand
`max(lithk - max(-topg, 0) rho_w/rho_i, 0)` on them against the model's
`max(h - h_f, 0)` on grounded cells, taken to the grid by the same operator.
Every pixel is weighed by its map-plane area, with af2 = 1 and without
maxmask1. The tool's own slvaf carries af2 and maxmask1 on top of this: that
is the area factor, cause 1 in `scalar_comparison_32km.md` (issue #97). The
+39.4 mm the 32 km run reached is the same map-plane term, cause 2 there, so
every number here compares with it like for like.
On the p2 run's 2015 state it reproduces the issue 13 attribution to the
tenth of a gigatonne: +1,594.9 Gt in total, -5,800.9 Gt of it in 4,455
grounding-line pixels.

The states are the v4 timing caches (BedMachine v4 cell averages on IU's
builds of each mesh, at 2015.0), the 32 km MAP, the 2 km Budd MAP of
24 September (snapshot 0241), the 32 km p2 and p4 annual output, and the
first year of the 2 km MAP-check control of 22 September. Each is changed two
ways, by the same rule on every mesh:

| change | what it measures |
|---|---|
| floating ice thinned by 25, 50 or 75 %, or removed | the tool's error alone: floating ice holds no VAF, so the model's VAF stays where it was |
| the shelves removed, and every marine grounded cell within X m of flotation with them | a grounding-line retreat onto open water; X is interpolated to the 550,000 km² of grounded area the 32 km run lost by 2300 |

## Result

Tool minus model slvaf in mm; positive is more sea-level rise in the tool's
reading. The mesh names the finest and the interior resolution. The timing
caches carry a 20 km ocean buffer, as the production layout does, and the
32 km meshes and the 2 km MAP's mesh carry none.

| mesh | cells | grounding-line pixels | shelves thinned 25 % | shelves removed | shelves removed and 550,000 km² retreat |
|---|---|---|---|---|---|
| 32 km, the p2 run in 2015 | 9,055 | 4,455 | +7.1 | +17.6 | +43.0 |
| 32 km MAP | 9,055 | 4,498 | +8.0 | +18.8 | +44.3 |
| 5 km, 50 km interior | 162,558 | 8,458 | +7.9 | +18.0 | +51.6 |
| 5 km, 100 km interior | 152,910 | 8,417 | +7.9 | +18.0 | +51.6 |
| 2.5 km, 25 km interior | 619,751 | 9,155 | +8.3 | +18.6 | +53.9 |
| 2.5 km, 50 km interior | 583,266 | 9,171 | +8.3 | +18.8 | +54.1 |
| 2 km, 20 km interior | 954,989 | 9,494 | +8.5 | +19.1 | +54.6 |
| 2 km, 40 km interior | 902,183 | 9,513 | +8.5 | +19.2 | +54.7 |
| 2 km MAP, 5 km interior | 1,835,718 | 7,016 | +8.2 | +17.9 | +54.1 |
| 1 km, 10 km interior (the production layout, IU's build) | 3,716,257 | 10,941 | +9.0 | +20.3 | +55.8 |
| 1 km, 20 km interior | 3,511,969 | 10,958 | +9.0 | +20.3 | +56.3 |
| 0.5 km, 10 km interior | 13,814,779 | 12,709 | +9.3 | +20.9 | +57.1 |
| 0.5 km, 5 km interior | 14,615,615 | 12,728 | +9.3 | +20.9 | +57.0 |

The first year of the 2 km control, on the MAP's mesh, gives the MAP's row
to 0.1 mm, so a forward's annual output and a MAP read alike.

The retreat stands in for the 32 km run's own evolution: on the run's 2015
state it gives +43.0 mm, and the run reached +39.4 mm. The run, p2 none:

| year | 2100 | 2150 | 2200 | 2250 | 2300 |
|---|---|---|---|---|---|
| model slvaf, mm | -186.2 | -356.9 | -467.2 | -468.3 | -373.6 |
| tool minus model, mm | +8.0 | +17.9 | +26.7 | +33.6 | +39.4 |
| floating ice, km³ | 517,500 | 325,900 | 149,600 | 57,200 | 16,100 |
| grounded area, million km² | 12.07 | 12.01 | 11.85 | 11.69 | 11.53 |

The error sits in the 8 km pixel means. A pixel that holds grounded ice
beside open water or thin shelf averages the two, and the flotation deficit
of the one cancels the VAF of the other. A finer mesh loses less in each
grounding-line pixel (with the shelves removed, 2.66 Gt at 32 km, 1.15 Gt at
5 km, 0.67 Gt at 1 km and 0.49 Gt at 0.5 km) and resolves more grounding
line, which crosses more pixels (4,455, 8,458, 10,941 and 12,709), and the
two cancel. The model's own VAF, from `limnsw`, carries none of it.

## For the decision

- At the production layout the tool's error for a given change in the ice
  sheet is the 32 km error or a little more: about +9 mm for shelves thinned
  by a quarter, +20 mm for shelves gone, +56 mm for the 32 km run's change
  by 2300.
- As a share, it depends on the production run's own slvaf. +56 mm is 1 to
  2 % of a contribution of 2.8 to 5.6 m, and 15 % of the 32 km run's slvaf
  at 2300, -373.6 mm.
- The model writes `limnsw` every year, so its own slvaf is at hand for the
  README beside the tool's. `grid_vaf_attribution.py` run on the production
  annual output gives the production run's own number.

## Provenance

The records and pixel planes are in IU's scratch directory
`ismip7_issue99/out` (purged after 30 days). Jobs, one core each, from the
checkout at `d1dced8` with the script before it was tracked: 10661220 and
10661244 (32 and 5 km), 10661224, 10661225 and 10661281 (2.5, 2 and 1 km),
10661248 (0.5 km). Job 10661353 ran the tracked script at `a3d52a6` through
`batch_runners/grid_vaf_attribution.script` on the 32 km p2 run in 2015 and
2300 and on the production layout, and reproduced every number. The
production layout loads in 3 min and builds its operator in 2.5 min, at
9 GB; a 0.5 km mesh takes 12 to 13 min and 6 to 7 min, at 24 GB.
