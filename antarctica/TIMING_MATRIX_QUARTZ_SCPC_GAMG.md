# Antarctica timing matrix

Generated: 2026-09-19 20:08 UTC

Campaign tag: `scpc_gamg_10step_dt0p125at2500_dg0_logvelnet_cached_strict_v4`. Configured lanes: 20; accepted: 16.

The primary timing covers only the 10-step transient loop; setup and checkpoint loading are recorded separately. The timestep is `0.125 × LC / 2500` years, capped at 0.125. Lanes use `scpc_gamg`, an exact-mesh prepared cache (prepared under `scpc_mumps`: the same initial state as the `scpc_mumps` campaign's lanes), and no rescue or subcycle recovery; the physics contract is `strict` (no apparent mass balance, no calving sink at the 2015 front). The initial state is the transferred 2.5 km MAP on every mesh (see Initial states).

## Time per year of simulation

Transient-loop wall time per simulated year, each lane at its own timestep (see Run status).

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (min/yr) | 32 cores (min/yr) | 64 cores (min/yr) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 32.3 | 15.7 | 10.0 |
| 1000 | 20000 | 1767017 | 3511969 | 25.1 | 14.2 | 8.9 |
| 2000 | 20000 | 483149 | 954989 | 3.8 | 2.0 | — |
| 2000 | 40000 | 456739 | 902183 | 3.8 | 1.9 | — |
| 2500 | 25000 | 314470 | 619751 | 2.0 | 1.1 | — |
| 2500 | 50000 | 296236 | 583266 | 1.9 | 1.0 | — |
| 5000 | 50000 | 83727 | 162558 | 0.5 | — | — |
| 5000 | 100000 | 78903 | 152910 | 0.5 | — | — |

### Time per 285-year simulation

The per-year cost extrapolated to 285 years, in hours (h) or days (d). It assumes the matrix timestep and the measured cost per step hold for the whole run, and leaves out setup, forcing updates and output.

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores | 32 cores | 64 cores |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 6.4 d | 3.1 d | 47.5 h |
| 1000 | 20000 | 1767017 | 3511969 | 5.0 d | 2.8 d | 42.3 h |
| 2000 | 20000 | 483149 | 954989 | 18.2 h | 9.3 h | — |
| 2000 | 40000 | 456739 | 902183 | 17.8 h | 8.9 h | — |
| 2500 | 25000 | 314470 | 619751 | 9.4 h | 5.0 h | — |
| 2500 | 50000 | 296236 | 583266 | 9.1 h | 4.7 h | — |
| 5000 | 50000 | 83727 | 162558 | 2.5 h | — | — |
| 5000 | 100000 | 78903 | 152910 | 2.4 h | — | — |

## Run status

| LC (m) | LC_coarse (m) | dt (yr) | 16 cores | 32 cores | 64 cores |
|---|---|---|---|---|---|
| 500 | 5000 | 0.025 | NOT PLANNED | NOT RUN | NOT RUN |
| 500 | 10000 | 0.025 | NOT PLANNED | NOT RUN | NOT RUN |
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
| 500 | 5000 | transferred 2.5 km MAP (policy) |
| 500 | 10000 | transferred 2.5 km MAP (policy) |
| 1000 | 10000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 1000 | 20000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 2000 | 20000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 2000 | 40000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 2500 | 25000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 2500 | 50000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 5000 | 50000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |
| 5000 | 100000 | campaign source MAP `inversion_icepack2_budd_n3_dg0_logvelnet_2500_25000_250iter.h5` |

### Status notes

- LC=500, LC_coarse=5000, 32 cores: **NOT RUN** — no record or status stamp.
- LC=500, LC_coarse=5000, 64 cores: **NOT RUN** — no record or status stamp.
- LC=500, LC_coarse=10000, 32 cores: **NOT RUN** — no record or status stamp.
- LC=500, LC_coarse=10000, 64 cores: **NOT RUN** — no record or status stamp.

## Successful timings

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (s) | 32 cores (s) | 64 cores (s) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 968.5 | 469.7 | 299.9 |
| 1000 | 20000 | 1767017 | 3511969 | 753.2 | 426.2 | 267.2 |
| 2000 | 20000 | 483149 | 954989 | 229.6 | 117.3 | — |
| 2000 | 40000 | 456739 | 902183 | 225.4 | 112.4 | — |
| 2500 | 25000 | 314470 | 619751 | 147.8 | 79.5 | — |
| 2500 | 50000 | 296236 | 583266 | 143.1 | 74.2 | — |
| 5000 | 50000 | 83727 | 162558 | 40.1 | — | — |
| 5000 | 100000 | 78903 | 152910 | 37.3 | — | — |

### Per-step timing

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (s/step) | 32 cores (s/step) | 64 cores (s/step) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 96.85 | 46.97 | 29.99 |
| 1000 | 20000 | 1767017 | 3511969 | 75.32 | 42.62 | 26.72 |
| 2000 | 20000 | 483149 | 954989 | 22.96 | 11.73 | — |
| 2000 | 40000 | 456739 | 902183 | 22.54 | 11.24 | — |
| 2500 | 25000 | 314470 | 619751 | 14.78 | 7.95 | — |
| 2500 | 50000 | 296236 | 583266 | 14.31 | 7.42 | — |
| 5000 | 50000 | 83727 | 162558 | 4.01 | — | — |
| 5000 | 100000 | 78903 | 152910 | 3.73 | — | — |

### Solver work

Newton iterations per step × iterations on the condensed velocity system per Newton iteration (back-substitutions under MUMPS, V-cycles under GAMG), averaged over the lane's diagnostic solves and counting the line search's solves. † marks an older record that holds only SNES's count, which leaves those solves out.

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores | 32 cores | 64 cores |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 12.9 × 34.3 | 12.9 × 34.1 | 12.9 × 34.2 |
| 1000 | 20000 | 1767017 | 3511969 | 11.5 × 35.2 | 11.5 × 38.4 | 11.5 × 38.2 |
| 2000 | 20000 | 483149 | 954989 | 10.9 × 36.5 | 10.9 × 35.5 | — |
| 2000 | 40000 | 456739 | 902183 | 11.2 × 33.5 | 11.2 × 34.0 | — |
| 2500 | 25000 | 314470 | 619751 | 10.7 × 40.7 | 10.7 × 39.5 | — |
| 2500 | 50000 | 296236 | 583266 | 12.4 × 32.0 | 12.4 × 32.6 | — |
| 5000 | 50000 | 83727 | 162558 | 11.4 × 43.7 | — | — |
| 5000 | 100000 | 78903 | 152910 | 10.1 × 47.4 | — | — |

Scout mass-residual acceptance limit: 5e-05 Gt.
