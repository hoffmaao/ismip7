r"""Run-shaping knobs: the discretization and physics choices that decide
which MAP a run writes, and which MAP a forward is allowed to load.

The inversion, the forward, the mesh pipeline, the probes, the preflight gate
and the per-core report all need these. Each used to re-declare its own
literal default, which is how ISMIP7_LC came to mean 2500 in ``simulation.py``
and ``preflight.py``, 8000 in ``inversion_icepack2.py`` and 32000 in
``thermo_prior.py``: with the variable unset, the gate blessed
``inversion_icepack2_budd_n3_dg0_2500.h5`` while the inversion would have
built and written the 8000 MAP. That is exactly the mismatch ``naming.py``
exists to prevent, so the values it builds names from are owned here and every
reader calls the accessors below rather than repeating a literal. A script
that genuinely needs a different value passes it at the call site or exports
the variable, so the divergence is visible instead of hiding in a default.

A knob left at its default is also absent from ``os.environ``, so
``core_report.py`` resolves the run-env block through this module: the report
is the only committed record of a run, so it has to state the value the run
used, not only the ones that happened to be exported.

Pure Python: importable without Firedrake so the preflight and the report stay
fast.
"""

import os

# 2500 m / 64 km is the production pair: it is the mesh the campaign inverts
# and runs on (``antarctica_64000_2500_buffered20000``), and the pair the
# README documents. The old 8000 and 32000 module-level defaults were
# dev-probe leftovers; a coarse probe now exports ISMIP7_LC / ISMIP7_LC_COARSE
# instead of disagreeing with the gate about what "unset" means.
LC_DEFAULT = "2500"
LC_COARSE_DEFAULT = "64000"
GEOMETRY_SPACE_DEFAULT = "dg0"
FRICTION_DEFAULT = "budd"
# THIS BRANCH (antarctica-n3) runs standard Glen n=3. An inversion and every
# forward that loads its MAP must agree on this.
N_FLOW_DEFAULT = "3.0"

GEOMETRY_SPACES = ("dg0", "cg1")


def lc():
    r"""Target edge length [m] in the refined region of the mesh."""
    return int(os.environ.get("ISMIP7_LC", LC_DEFAULT))


def lc_coarse():
    r"""Target edge length [m] in the coarse region of the mesh."""
    return int(os.environ.get("ISMIP7_LC_COARSE", LC_COARSE_DEFAULT))


def geometry_space():
    r"""Discretization of h/s/b, ``'dg0'`` or ``'cg1'``. Validated here so
    every reader rejects the same set."""
    value = os.environ.get(
        "ISMIP7_GEOMETRY_SPACE", GEOMETRY_SPACE_DEFAULT).lower()
    if value not in GEOMETRY_SPACES:
        raise ValueError(
            f"ISMIP7_GEOMETRY_SPACE must be 'dg0' or 'cg1', got {value!r}"
        )
    return value


def friction():
    r"""Friction law: ``budd``, ``regularized_coulomb`` or ``budd_legacy``."""
    return os.environ.get("ISMIP7_FRICTION", FRICTION_DEFAULT)


def n_flow():
    r"""Glen flow-law exponent."""
    return float(os.environ.get("ISMIP7_N_FLOW", N_FLOW_DEFAULT))


# Calving front (icepack2_tools.levelset). ``none`` is the pre-Sep-2026
# behaviour: on a buffered mesh the front advances freely and never calves.
CALVING_DEFAULT = "none"
CALVING_LAWS = ("none", "fixed", "vonmises")
# ISSM defaults for the von Mises thresholds (Morlighem et al. 2016).
CALVING_SIGMA_MAX_GROUNDED_DEFAULT = "1.0"     # MPa
CALVING_SIGMA_MAX_FLOATING_DEFAULT = "0.15"    # MPa


def calving_law():
    r"""``ISMIP7_CALVING``: ``none``, ``fixed`` or ``vonmises``."""
    value = os.environ.get("ISMIP7_CALVING", CALVING_DEFAULT).lower()
    if value not in CALVING_LAWS:
        raise ValueError(
            f"ISMIP7_CALVING must be one of {CALVING_LAWS}, got {value!r}"
        )
    return value


def fixed_front():
    r"""``ISMIP7_FIXED_FRONT``: the legacy pinned front.

    On when the variable is set to anything but the exact string ``"0"``:
    ``run_core_matrix.sh`` exports it unconditionally, so ``=0`` has to be the
    way to turn it off from there.
    """
    return os.environ.get("ISMIP7_FIXED_FRONT") not in (None, "0")


def calving_sigma_max():
    r"""Von Mises thresholds (grounded, floating) [MPa]."""
    return (
        float(os.environ.get("ISMIP7_CALVING_SIGMA_MAX_GROUNDED",
                             CALVING_SIGMA_MAX_GROUNDED_DEFAULT)),
        float(os.environ.get("ISMIP7_CALVING_SIGMA_MAX_FLOATING",
                             CALVING_SIGMA_MAX_FLOATING_DEFAULT)),
    )
