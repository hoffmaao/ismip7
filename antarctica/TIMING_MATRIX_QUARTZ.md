# Antarctica timing matrix

Generated: 2026-09-18 14:06 UTC

Campaign tag: `scpc_mumps_10step_dt0p125at2500_dg0_logvelnet_cached_strict_v4`. Configured lanes: 20; accepted: 16.

The primary timing covers only the 10-step transient loop; setup and checkpoint loading are recorded separately. The timestep is `0.125 × LC / 2500` years, capped at 0.125. Lanes use `scpc_mumps`, an exact-mesh prepared cache, and no rescue or subcycle recovery; the physics contract is `strict` (no apparent mass balance, no calving sink at the 2015 front). The initial state is the transferred 2.5 km MAP on every mesh (see Initial states).

## Time per year of simulation

Transient-loop wall time per simulated year, each lane at its own timestep (see Run status).

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (min/yr) | 32 cores (min/yr) | 64 cores (min/yr) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 107.9 | 72.5 | 56.7 |
| 1000 | 20000 | 1767017 | 3511969 | 76.0 | 46.6 | 37.3 |
| 2000 | 20000 | 483149 | 954989 | 10.4 | 6.3 | — |
| 2000 | 40000 | 456739 | 902183 | 9.7 | 7.3 | — |
| 2500 | 25000 | 314470 | 619751 | 4.4 | 2.8 | — |
| 2500 | 50000 | 296236 | 583266 | 5.8 | 3.7 | — |
| 5000 | 50000 | 83727 | 162558 | 1.1 | — | — |
| 5000 | 100000 | 78903 | 152910 | 0.9 | — | — |

### Time per 285-year simulation

The per-year cost extrapolated to 285 years, in hours (h) or days (d). It assumes the matrix timestep and the measured cost per step hold for the whole run, and leaves out setup, forcing updates and output.

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores | 32 cores | 64 cores |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 21.3 d | 14.4 d | 11.2 d |
| 1000 | 20000 | 1767017 | 3511969 | 15.1 d | 9.2 d | 7.4 d |
| 2000 | 20000 | 483149 | 954989 | 2.1 d | 29.8 h | — |
| 2000 | 40000 | 456739 | 902183 | 46.2 h | 34.8 h | — |
| 2500 | 25000 | 314470 | 619751 | 20.8 h | 13.5 h | — |
| 2500 | 50000 | 296236 | 583266 | 27.5 h | 17.4 h | — |
| 5000 | 50000 | 83727 | 162558 | 5.3 h | — | — |
| 5000 | 100000 | 78903 | 152910 | 4.3 h | — | — |

## Run status

| LC (m) | LC_coarse (m) | dt (yr) | 16 cores | 32 cores | 64 cores |
|---|---|---|---|---|---|
| 500 | 5000 | 0.025 | NOT PLANNED | NOT RUN | NOT RUN |
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
| 500 | 5000 | transferred 2.5 km MAP (policy) |
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

- LC=500, LC_coarse=5000, 32 cores: **NOT RUN** — no record or status stamp.
- LC=500, LC_coarse=5000, 64 cores: **NOT RUN** — no record or status stamp.
- LC=500, LC_coarse=10000, 32 cores: **NUMERICAL FAILURE** — runaway_tripwire at step-1.
- LC=500, LC_coarse=10000, 64 cores: **BLOCKED BY SCOUT** — scout did not pass.

## Successful timings

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (s) | 32 cores (s) | 64 cores (s) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 3235.9 | 2176.5 | 1701.1 |
| 1000 | 20000 | 1767017 | 3511969 | 2281.3 | 1397.1 | 1119.2 |
| 2000 | 20000 | 483149 | 954989 | 626.2 | 376.1 | — |
| 2000 | 40000 | 456739 | 902183 | 584.0 | 440.0 | — |
| 2500 | 25000 | 314470 | 619751 | 328.3 | 213.2 | — |
| 2500 | 50000 | 296236 | 583266 | 433.6 | 275.3 | — |
| 5000 | 50000 | 83727 | 162558 | 84.0 | — | — |
| 5000 | 100000 | 78903 | 152910 | 67.8 | — | — |

### Per-step timing

| LC (m) | LC_coarse (m) | Vertices | Cells | 16 cores (s/step) | 32 cores (s/step) | 64 cores (s/step) |
|---|---|---|---|---|---|---|
| 500 | 5000 | — | — | — | — | — |
| 500 | 10000 | — | — | — | — | — |
| 1000 | 10000 | 1869088 | 3716257 | 323.59 | 217.65 | 170.11 |
| 1000 | 20000 | 1767017 | 3511969 | 228.13 | 139.71 | 111.92 |
| 2000 | 20000 | 483149 | 954989 | 62.62 | 37.61 | — |
| 2000 | 40000 | 456739 | 902183 | 58.40 | 44.00 | — |
| 2500 | 25000 | 314470 | 619751 | 32.83 | 21.32 | — |
| 2500 | 50000 | 296236 | 583266 | 43.36 | 27.53 | — |
| 5000 | 50000 | 83727 | 162558 | 8.40 | — | — |
| 5000 | 100000 | 78903 | 152910 | 6.78 | — | — |

Scout mass-residual acceptance limit: 5e-05 Gt.
