# Antarctica timing matrix

Generated: 2026-09-19 19:33 UTC

Campaign tag: `scpc_mumps_10step_dt0p125at2500_dg0_logvelnet_cached_strict_v4`. Configured lanes: 20; accepted: 16.

The primary timing covers only the 10-step transient loop; setup and checkpoint loading are recorded separately. The timestep is `0.125 × LC / 2500` years, capped at 0.125. Lanes use `scpc_mumps`, an exact-mesh prepared cache, and no rescue or subcycle recovery; the physics contract is `strict` (no apparent mass balance, no calving sink at the 2015 front). The initial state is the transferred 2.5 km MAP on every mesh (see Initial states).

## Time per year of simulation

Transient-loop wall time per simulated year, each lane at its own timestep (see Run status).

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (min/yr) | 32 cores (min/yr) | 64 cores (min/yr) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 34.0 | 23.4 | 23.2 |
| 1000 | 20000 | 1767017 | 3511969 | 25.8 | 16.8 | 15.2 |
| 2000 | 20000 | 483149 | 954989 | 3.2 | 2.0 | — |
| 2000 | 40000 | 456739 | 902183 | 2.6 | 2.0 | — |
| 2500 | 25000 | 314470 | 619751 | 1.6 | 1.1 | — |
| 2500 | 50000 | 296236 | 583266 | 1.6 | 1.2 | — |
| 5000 | 50000 | 83727 | 162558 | 0.4 | — | — |
| 5000 | 100000 | 78903 | 152910 | 0.3 | — | — |

### Time per 285-year simulation

The per-year cost extrapolated to 285 years, in hours (h) or days (d). It assumes the matrix timestep and the measured cost per step hold for the whole run, and leaves out setup, forcing updates and output.

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores | 32 cores | 64 cores |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 6.7 d | 4.6 d | 4.6 d |
| 1000 | 20000 | 1767017 | 3511969 | 5.1 d | 3.3 d | 3.0 d |
| 2000 | 20000 | 483149 | 954989 | 15.1 h | 9.3 h | — |
| 2000 | 40000 | 456739 | 902183 | 12.2 h | 9.4 h | — |
| 2500 | 25000 | 314470 | 619751 | 7.4 h | 5.1 h | — |
| 2500 | 50000 | 296236 | 583266 | 7.7 h | 5.5 h | — |
| 5000 | 50000 | 83727 | 162558 | 2.0 h | — | — |
| 5000 | 100000 | 78903 | 152910 | 1.6 h | — | — |

## Run status

| LC (m) | LC_coarse (m) | dt (yr) | 16 cores | 32 cores | 64 cores |
|---|---|---|---|---|---|
| 500 | 5000 | 0.025 | NOT PLANNED | NUMERICAL FAILURE | NOT RUN |
| 500 | 10000 | 0.025 | NOT PLANNED | NUMERICAL FAILURE | BLOCKED BY SCOUT |
| 1000 | 10000 | 0.05 | OK | OK | OK |
| 1000 | 20000 | 0.05 | OK | OK | OK |
| 2000 | 20000 | 0.1 | OK | OK | NOT PLANNED |
| 2000 | 40000 | 0.1 | OK | OK | NOT PLANNED |
| 2500 | 25000 | 0.125 | OK | OK | NOT PLANNED |
| 2500 | 50000 | 0.125 | OK | OK | NOT PLANNED |
| 5000 | 50000 | 0.125 | OK | NOT PLANNED | NOT PLANNED |
| 5000 | 100000 | 0.125 | OK | NOT PLANNED | NOT PLANNED |

### Initial states

| LC (m) | LC_coarse (m) | Initial state |
|---|---|---|
| 500 | 5000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 500 | 10000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 1000 | 10000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 1000 | 20000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 2000 | 20000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 2000 | 40000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 2500 | 25000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 2500 | 50000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 5000 | 50000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 5000 | 100000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |

### Status notes

- LC=500, LC_coarse=5000, 32 cores: **NUMERICAL FAILURE** — runaway_tripwire at step-1.
- LC=500, LC_coarse=5000, 64 cores: **NOT RUN** — no record or status stamp.
- LC=500, LC_coarse=10000, 32 cores: **NUMERICAL FAILURE** — runaway_tripwire at step-1.
- LC=500, LC_coarse=10000, 64 cores: **BLOCKED BY SCOUT** — scout did not pass.

## Successful timings

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (s) | 32 cores (s) | 64 cores (s) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 1019.7 | 702.5 | 694.9 |
| 1000 | 20000 | 1767017 | 3511969 | 774.9 | 505.0 | 456.0 |
| 2000 | 20000 | 483149 | 954989 | 191.3 | 117.3 | — |
| 2000 | 40000 | 456739 | 902183 | 154.2 | 118.7 | — |
| 2500 | 25000 | 314470 | 619751 | 116.6 | 80.3 | — |
| 2500 | 50000 | 296236 | 583266 | 121.2 | 86.4 | — |
| 5000 | 50000 | 83727 | 162558 | 31.5 | — | — |
| 5000 | 100000 | 78903 | 152910 | 24.9 | — | — |

### Per-step timing

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (s/step) | 32 cores (s/step) | 64 cores (s/step) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 101.97 | 70.25 | 69.49 |
| 1000 | 20000 | 1767017 | 3511969 | 77.49 | 50.50 | 45.60 |
| 2000 | 20000 | 483149 | 954989 | 19.13 | 11.73 | — |
| 2000 | 40000 | 456739 | 902183 | 15.42 | 11.87 | — |
| 2500 | 25000 | 314470 | 619751 | 11.66 | 8.03 | — |
| 2500 | 50000 | 296236 | 583266 | 12.12 | 8.64 | — |
| 5000 | 50000 | 83727 | 162558 | 3.15 | — | — |
| 5000 | 100000 | 78903 | 152910 | 2.49 | — | — |

### Solver work

Newton iterations per step × iterations on the condensed velocity system per Newton iteration (back-substitutions under MUMPS, V-cycles under GAMG), averaged over the lane's diagnostic solves and counting the line search's solves. † marks an older record that holds only SNES's count, which leaves those solves out.

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores | 32 cores | 64 cores |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 12.9 × 2.2 | 12.9 × 2.2 | 12.9 × 2.2 |
| 1000 | 20000 | 1767017 | 3511969 | 11.4 × 2.1 | 11.4 × 2.1 | 11.4 × 2.1 |
| 2000 | 20000 | 483149 | 954989 | 10.9 × 2.1 | 10.9 × 2.1 | — |
| 2000 | 40000 | 456739 | 902183 | 11.2 × 2.2 | 11.2 × 2.2 | — |
| 2500 | 25000 | 314470 | 619751 | 10.7 × 2.1 | 10.7 × 2.1 | — |
| 2500 | 50000 | 296236 | 583266 | 12.4 × 2.2 | 12.4 × 2.2 | — |
| 5000 | 50000 | 83727 | 162558 | 11.4 × 2.1 | — | — |
| 5000 | 100000 | 78903 | 152910 | 10.1 × 2.2 | — | — |

Scout mass-residual acceptance limit: 5e-05 Gt.
