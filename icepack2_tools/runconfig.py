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

# 1000 m / 10 km is the production pair (``antarctica_10000_1000_buffered20000``):
# the finest mesh the Quartz timing matrix carries through a 285-year run in
# two days, under scpc_gamg on 64 ranks
# (antarctica/TIMING_MATRIX_QUARTZ_SCPC_GAMG.md). It is the pair the batch
# runners export (batch_runners/site_env.sh) and the README documents; until
# 2026-09-19 that was 2500 m / 64 km. The old 8000 and 32000 module-level
# defaults were dev-probe leftovers; a coarse probe exports ISMIP7_LC /
# ISMIP7_LC_COARSE instead of disagreeing with the gate about what "unset" means.
LC_DEFAULT = "1000"
LC_COARSE_DEFAULT = "10000"
GEOMETRY_SPACE_DEFAULT = "dg0"
FRICTION_DEFAULT = "budd"
# THIS BRANCH (antarctica-n3) runs standard Glen n=3. An inversion and every
# forward that loads its MAP must agree on this.
N_FLOW_DEFAULT = "3.0"
# How the Budd law (dual_friction.budd_nhat) zeroes shelf friction. Before
# 2026-09-13 it tested the sign of N = max(p_I - p_W, 0), a roundoff residue
# on floating ice, and the delta floor lifted every roundoff-positive shelf
# cell to the friction cap (hoffmaao/antarctica e602705). "haf" gates on
# height above flotation. Stamped into every Budd state checkpoint and required
# by the timing-cache manifest, so a state solved under the old gate can never
# seed a fixed-law lane.
BUDD_SHELF_GATE = "haf"

# Exact-mesh timing caches initialize DG0 geometry from BedMachine on the
# TARGET mesh. The raster is first sampled into CG1 and then L2-projected to
# DG0, so each cell stores an average rather than one centroid pixel. Keep the
# method name stable: it is stamped into cache provenance and changing the
# construction must invalidate old caches.
TARGET_MESH_GEOMETRY_METHOD = "target-native-bedmachine-cell-average-v1"

GEOMETRY_SPACES = ("dg0", "cg1")

# How a raster (BedMachine) is put onto a DG0 geometry cell.
#   vertex    - icepack's bilinear sample at the three CG1 vertices, then the
#               L2 projection of that linear interpolant (= the mean of the
#               3 vertex values). The pre-Sep-2026 behaviour. A 20 km interior
#               cell sees 3 of its ~1600 BedMachine pixels.
#   cell_mean - the mean of the raster over the cell itself, sampled on an
#               equal-area sub-triangle lattice at pixel density
#               (geometry.raster_cell_mean).
# MAPs record the method used; the forward reads it back from the MAP.
RASTER_SAMPLES = ("vertex", "cell_mean")
RASTER_SAMPLE_DEFAULT = "vertex"


# Floor-cell coercivity drags of the dual-friction residual (dual_friction.py
# ``ocean_drag`` / ``u_lim``). The forward has applied them since Jul 2026; the
# inversion never passed them, so a MAP that solved the inversion's F was not
# a solution of the forward's F in the ice-free buffer cells: the same mixed
# state measured ||F||=1.3e1 in the inversion and 1.5e10 in the forward
# (2026-09-14), and every strict scout that trusted it ran away within four
# steps. Both now read the knobs here. ``u_lim`` is only the threshold of the
# soft speed limiter; its gain ``k_lim`` stays a live Constant at 0 in the
# forward (raised for rescue solves only) and is passed as 0 by the inversion.
OCEAN_DRAG_DEFAULT = "1e-2"   # MPa yr/m, linear drag ramping to zero at H_OCEAN
H_OCEAN_DEFAULT = "10.0"      # m
U_LIM_DEFAULT = "2e4"         # m/yr, ~5x the fastest observed Antarctic flow


def residual_stabilizers():
    r"""``ocean_drag``, ``h_ocean`` and ``u_lim`` keyword arguments for
    ``build_rc_residual``; identical in the inversion and the forward."""
    return {
        "ocean_drag": float(os.environ.get("ISMIP7_OCEAN_DRAG", OCEAN_DRAG_DEFAULT)),
        "h_ocean": float(os.environ.get("ISMIP7_H_OCEAN", H_OCEAN_DEFAULT)),
        "u_lim": float(os.environ.get("ISMIP7_U_LIM", U_LIM_DEFAULT)),
    }


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


def raster_sample():
    r"""How BedMachine is sampled onto a DG0 cell: ``'vertex'`` or
    ``'cell_mean'``. See RASTER_SAMPLES."""
    value = os.environ.get(
        "ISMIP7_RASTER_SAMPLE", RASTER_SAMPLE_DEFAULT).lower()
    if value not in RASTER_SAMPLES:
        raise ValueError(
            f"ISMIP7_RASTER_SAMPLE must be one of {RASTER_SAMPLES}, got {value!r}"
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


FRACTURE_MODES = ("none", "mask")
FRACTURE_DEFAULT = "none"


def fracture():
    r"""``ISMIP7_FRACTURE``: how the ISMIP7 ice-shelf collapse forcing is
    applied. ``none`` (default) loads nothing; ``mask`` removes the ice of
    every FLOATING cell the year's collapse mask flags, booked as calving
    (protocol path C, discussions #30 and #33: floating ice only; no mask
    exists for historical or OCX, so those runs see nothing either way).
    A stress-gated variant (Lai et al. 2020) is not implemented."""
    value = os.environ.get("ISMIP7_FRACTURE", FRACTURE_DEFAULT).lower()
    if value not in FRACTURE_MODES:
        raise ValueError(f"ISMIP7_FRACTURE must be one of {FRACTURE_MODES}, got {value!r}")
    return value


def ismip7_output():
    r"""``ISMIP7_OUTPUT``: record the ISMIP7 yearly fields and scalars.

    ``1`` enables it; ``0`` and the empty string disable it. The value set is
    closed, like ``fracture`` and ``apparent_mb_mode``: there is one spelling
    each way, and anything else raises rather than silently deciding whether
    a submission gets written.
    """
    value = (os.environ.get("ISMIP7_OUTPUT") or "").strip()
    if value in ("", "0"):
        return False
    if value == "1":
        return True
    raise ValueError(
        f"ISMIP7_OUTPUT must be 1 to enable or 0/empty to disable, "
        f"got {value!r}"
    )


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


def apparent_mb_mode():
    r"""``ISMIP7_APPARENT_MB``: the apparent-mass-balance init, or None for off.

    ``"div"`` cancels only the flux divergence, so the t=0 tendency is
    SMB minus melt (gia-style). ``"1"`` or ``"balance"`` also subtracts the
    initial forcing, so the t=0 tendency is exactly zero: a balanced control
    in the ISMIP6 ctrl_proj sense.

    ``0``, ``off``, ``none`` and the empty string mean OFF. The batch runners
    export this unconditionally and ``sbatch --export=ALL,VAR=...`` cannot
    unset a variable, so there has to be an off value; without one a run asked
    to drop the correction would silently get the full balanced one.

    The value set is closed, like ``calving_law`` and ``raster_sample``:
    anything else raises. ``no`` and ``false`` are not off spellings, and
    ``divergence`` is not ``div``, so accepting them would hand back the
    balanced control, which differs from both by the whole t=0 forcing.
    """
    value = (os.environ.get("ISMIP7_APPARENT_MB") or "").strip().lower()
    if value in ("", "0", "off", "none"):
        return None
    if value in ("1", "balance"):
        return "balance"
    if value == "div":
        return "div"
    raise ValueError(
        f"ISMIP7_APPARENT_MB must be 1 or balance (balanced control), div "
        f"(divergence only), or 0/off/none/empty to disable; got {value!r}"
    )


def auto_resume():
    r"""``ISMIP7_AUTO_RESUME``: continue unattended from this experiment's own
    newest checkpoint when no explicit restart is given.

    An integer flag, so ``=0`` turns it OFF. The batch runners export it
    unconditionally and ``sbatch --export=ALL,VAR=...`` gives no way to unset a
    variable, so ``0`` has to be the off switch; testing the string for mere
    presence would silently resume a run the user asked to start clean.
    """
    value = (os.environ.get("ISMIP7_AUTO_RESUME") or "").strip()
    if not value:
        return False
    try:
        return int(value) != 0
    except ValueError:
        raise ValueError(
            f"ISMIP7_AUTO_RESUME must be an integer flag (0 to disable), "
            f"got {value!r}"
        ) from None


def calving_sigma_max():
    r"""Von Mises thresholds (grounded, floating) [MPa]."""
    return (
        float(os.environ.get("ISMIP7_CALVING_SIGMA_MAX_GROUNDED",
                             CALVING_SIGMA_MAX_GROUNDED_DEFAULT)),
        float(os.environ.get("ISMIP7_CALVING_SIGMA_MAX_FLOATING",
                             CALVING_SIGMA_MAX_FLOATING_DEFAULT)),
    )


# ── Roots a second checkout does not carry ──────────────────────────────
# The code moves with the invocation; the large gitignored artifacts do not.
# site_env.sh names them for the shell half of a run; these are the Python
# half, so a driver started by hand resolves the same paths.
_ANTARCTICA = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "antarctica"
)


def obs_data_root():
    r"""``ISMIP7_OBS_DATA_ROOT``: BedMachine, MEaSUREs, RACMO and the dH/dt
    cache rasters. Falls back to this checkout's ``antarctica/data``, which is
    right wherever the data sits beside the code and empty where it does not."""
    return os.environ.get("ISMIP7_OBS_DATA_ROOT",
                          os.path.join(_ANTARCTICA, "data"))


def k_per_basin_candidates(results_dir, lc_value):
    r"""Where to look for the calibrated per-basin K, in order.

    ``ISMIP7_K_PER_BASIN_NPZ`` wins; then this mesh's calibration and the
    2500 m fallback (16 basin scalars remapped through the IMBIE2 8 km grid, so
    mesh-independent). The file is gitignored, so a checkout that does not carry
    it warns and falls back to an SMB-only melt source; name it with the
    override when it lives elsewhere.
    """
    override = os.environ.get("ISMIP7_K_PER_BASIN_NPZ")
    if override:
        return [override]
    names = [f"calibrated_K_per_basin_{lc_value}.npz",
             "calibrated_K_per_basin_2500.npz"]
    # At lc=2500 the two names coincide.
    seen, out = set(), []
    for name in names:
        path = os.path.join(results_dir, name)
        if path not in seen:
            seen.add(path)
            out.append(path)
    return out
