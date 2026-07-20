# Branch `antarctica-n3`: standard Glen n=3

This branch is a parallel line to `antarctica`, identical in every respect
**except the ice flow-law exponent**: it runs the composite viscous rheology
at **n = 3 (standard Glen)** instead of n = 4 (Goldsby-Kohlstedt dislocation
creep). n = 3 is the ISMIP6/ISMIP7 default and the more directly comparable
configuration; n = 4 remains on `antarctica` untouched.

## What changed vs `antarctica`

The n = 4 assumption lived entirely in two environment-variable defaults, so
the branch difference is small and self-documenting:

| knob | `antarctica` (n=4) | `antarctica-n3` (n=3) | why |
|---|---|---|---|
| `ISMIP7_N_FLOW` | `4.0` | `3.0` | the flow exponent |
| `ISMIP7_A4_FACTOR` | `10.0` | `1.0` | `A0 = rate_factor(260 K)` **is** the n=3 fluidity, so the composite main term is plain Glen and needs no rescale. At n=4 the factor lifts `A_3` to `A_4` so `A_4·τc⁴ ≈ A_3·τc³` at `τc`. |

Both defaults are module constants (`N_FLOW_DEFAULT`, `A4_FACTOR_DEFAULT`) in
`simulation.py` and `inversion_icepack2.py`; a run still overrides them with
the env vars. The inversion and the forward that loads its MAP **must** use
the same pair.

## MAPs

The inversion is n-specific (the fluidity, the `τc^{n-1}` linearization, and
the n-continuation target all depend on n), so n = 3 needs its **own** MAPs.
They carry an `_n3` filename tag so they coexist on disk with the n = 4 MAPs:

- n = 4: `inversion_icepack2_budd_<lc>.h5` (untagged, legacy)
- n = 3: `inversion_icepack2_budd_n3_<lc>.h5`

The tag is produced by `map_n_tag()` (empty at n = 4 for backward
compatibility). MAP h5 files are gitignored and regenerated per machine:

```
# 32 km n=3 Budd MAP (gated launcher; logs -> results/logs/budd_inv_n3_*)
ISMIP7_FRICTION=budd ISMIP7_LC=32000 ISMIP7_N_FLOW=3.0 \
  ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_320000_32000.msh \
  NRANKS=16 MIN_FREE_GB=64 MIN_FREE_CORES=20 \
  setsid nohup antarctica/scripts/launch_rc_inversion_when_ready.sh &
```

Once the MAP exists, the forward path (CTRL, projections, the core-experiment
matrix) runs exactly as on `antarctica` - `ISMIP7_N_FLOW`/`ISMIP7_A4_FACTOR`
carry through `setup_model`, and the launcher/log/lock names are n-tagged so
n = 3 and n = 4 inversions can run on the same machine.

See `COMPOSITE_RHEOLOGY.md` for the full formulation.
