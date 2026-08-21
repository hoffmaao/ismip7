r"""Reference-climate window and pool scenario for the ISMIP7 runs.

The control's SMB climatology and the projections' aSMB re-reference pool must
be built from the SAME (window, scenario) pair. If they disagree, projection
minus control differences two different baselines - the same class of error as
branching a control from a cold start instead of the historical.

The control run, the projection driver, the preflight gate and the per-core
report all need these values, so they are owned here instead of being copied
into each reader. A copy that drifts silently reintroduces the baseline
mismatch, and the report - the only committed record of a run - would then
document an effective value the run never used.

The pool scenario is ssp126, per the ISMIP7 protocol cheat sheet (April 2026):
"The climatology (2000-2029) should be created using a combination of the
historical and the SSP126 simulations" (discussion #28: last 15 yr historical
+ first 15 yr ssp126).

Pure Python: importable without Firedrake, so the preflight and the report
stay fast.
"""

import os

CLIM_START_DEFAULT = "2000"
CLIM_END_DEFAULT = "2029"
CLIM_SCENARIO_DEFAULT = "ssp126"


def clim_start():
    r"""First year of the reference-climate window."""
    return int(os.environ.get("ISMIP7_CLIM_START", CLIM_START_DEFAULT))


def clim_end():
    r"""Last year of the reference-climate window (inclusive)."""
    return int(os.environ.get("ISMIP7_CLIM_END", CLIM_END_DEFAULT))


def clim_scenario():
    r"""Scenario pooled with ``historical`` to build the climatology."""
    return os.environ.get("ISMIP7_CLIM_SCENARIO", CLIM_SCENARIO_DEFAULT)
