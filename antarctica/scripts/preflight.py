#!/usr/bin/env python3
r"""Data preflight for the ISMIP7 core experiments.

Answers "which experiments can run on this machine right now?" in a few
seconds, checking every input each core needs: mesh, MAP inversion,
boundary ids, the melt calibration, RACMO, the OI climatology (core 11's stopgap),
and the (ESM, scenario) atmosphere/ocean trees over the run period, the
control's `ctrl` ocean among them. Honors the same
environment knobs as the runs (ISMIP7_LC, ISMIP7_FRICTION,
ISMIP7_OI_VERSION, ...). The run period is read from each core's driver,
with ISMIP7_T_START and ISMIP7_T_END applied as that driver applies them, and
a cold start that the geometry backdating refuses (issue 117) is BLOCKED.

Usage:
    python scripts/preflight.py
    ISMIP7_LC=500 ISMIP7_FRICTION=regularized_coulomb python scripts/preflight.py
"""
import ast, glob, math, os, sys

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_ANT = os.path.dirname(_SCRIPTS)
_PROJECT = os.path.dirname(_ANT)
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from icepack2_tools.forcing import (
    ISMIP7Atmosphere, ISMIP7Ocean, OCX, OCX_ATMOSPHERE_SOURCE,
    _oi_climatology_path, _find_ismip7_data, imbie2_basin_path, forcing_year,
)
from icepack2_tools.boundary import sidecar_path
from icepack2_tools.naming import map_basename
from icepack2_tools.climatology import clim_start, clim_end, clim_scenario
from icepack2_tools.runconfig import (
    obs_data_root,
    calving_law as _calving_law, calving_sigma_max as _calving_sigma_max,
    friction as _friction, geometry_space as _geometry_space, lc as _lc,
    lc_coarse as _lc_coarse, ocx_forcing as _ocx_forcing, ocx_ocean as _ocx_ocean,
    mesh_override, deltat_per_basin_npz, k_per_basin_npz,
    melt_calibration_contract, geometry_backdate_years,
)
DATA_DIR = obs_data_root()
from mesh_naming import get_buffer_m, mesh_filename

MESH_DIR = os.path.join(_ANT, "mesh")

lc = _lc()
lc_coarse = _lc_coarse()
friction = _friction()
oi_version = os.environ.get("ISMIP7_OI_VERSION", "30_sep")
CLIM_START = clim_start()
CLIM_END = clim_end()
CLIM_SCENARIO = clim_scenario()
root = _find_ismip7_data()

#: core, title, the driver under scripts/ that runs it, ESM, scenario. The
#: years a core is checked over come from its driver (driver_period).
CORES = [
    (1, "hist CESM2-WACCM", "historical/cesm_waccm.py", "CESM2-WACCM", "historical"),
    (2, "hist MRI-ESM2-0", "historical/mri_esm2.py", "MRI-ESM2-0", "historical"),
    (3, "ssp370 CESM2-WACCM", "projections/ssp370_cesm_waccm.py", "CESM2-WACCM", "ssp370"),
    (4, "ssp370 MRI-ESM2-0", "projections/ssp370_mri_esm2.py", "MRI-ESM2-0", "ssp370"),
    (5, "ssp126 CESM2-WACCM", "projections/ssp126_cesm_waccm.py", "CESM2-WACCM", "ssp126"),
    (6, "ssp126 MRI-ESM2-0", "projections/ssp126_mri_esm2.py", "MRI-ESM2-0", "ssp126"),
    (7, "ssp585 CESM2-WACCM", "projections/ssp585_cesm_waccm.py", "CESM2-WACCM", "ssp585"),
    (8, "ssp585 MRI-ESM2-0", "projections/ssp585_mri_esm2.py", "MRI-ESM2-0", "ssp585"),
    (9, "CTRL2015 (CESM2-WACCM clim)", "control/run.py", "CESM2-WACCM", None),
    (10, "CTRL2015 (MRI-ESM2-0 clim)", "control/run.py", "MRI-ESM2-0", None),
    (11, "OCX obs-constrained", "projections/ocx.py", None, None),
]


def _keyword(call, name, driver):
    r"""The literal a driver passes ``run_core_experiment`` as ``name``."""
    for k in call.keywords:
        if k.arg == name:
            return ast.literal_eval(k.value)
    raise ValueError(f"{driver} passes run_core_experiment no {name}")


def _time(node, env):
    r"""A start or end time as a driver computes it: a number, or
    ``float(os.environ.get(knob, default))`` evaluated under ``env``."""
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if (isinstance(node, ast.Call) and ast.unparse(node.func) == "float"
            and len(node.args) == 1 and isinstance(node.args[0], ast.Call)
            and ast.unparse(node.args[0].func) == "os.environ.get"
            and len(node.args[0].args) == 2):
        knob, default = (ast.literal_eval(a) for a in node.args[0].args)
        return float(env.get(knob, default))
    raise ValueError(f"cannot read a time from {ast.unparse(node)}")


def driver_period(driver, core=None, env=None):
    r"""``(t_start, t_end)`` of the run ``driver`` makes under ``env``
    (``os.environ`` by default), read from its source: importing a driver
    loads Firedrake, which takes about 20 s.

    A core-experiment shim passes ``t_start_default`` and ``t_end_default``
    to ``experiment.run_core_experiment``, where ``ISMIP7_T_START`` and
    ``ISMIP7_T_END`` replace them. OCX and the control set module-level
    ``T_START`` and ``T_END``, each a number or ``float(os.environ.get(knob,
    default))``. Any other form raises ``ValueError``, as does a shim that
    runs a core other than ``core``, and ``main`` reports that core BLOCKED
    with the reason.
    """
    env = os.environ if env is None else env
    path = os.path.join(_SCRIPTS, driver)
    with open(path) as f:
        tree = ast.parse(f.read(), path)
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and ast.unparse(node.func) == "run_core_experiment":
            runs = _keyword(node, "core", driver)
            if core is not None and runs != core:
                raise ValueError(f"{driver} runs core {runs}, listed here as core {core}")
            return (float(env.get("ISMIP7_T_START", _keyword(node, "t_start_default", driver))),
                    float(env.get("ISMIP7_T_END", _keyword(node, "t_end_default", driver))))
    times = {target.id: _time(node.value, env)
             for node in tree.body if isinstance(node, ast.Assign)
             for target in node.targets
             if isinstance(target, ast.Name) and target.id in ("T_START", "T_END")}
    if set(times) != {"T_START", "T_END"}:
        raise ValueError(f"{driver} has no run_core_experiment call and no "
                         f"module-level T_START and T_END")
    return times["T_START"], times["T_END"]


def atm_years(esm, scenario, var="acabf-anomaly"):
    r"""The years the reader can open, by the reader's own match: a gate that
    counts files the loader would not find clears a run that then fails."""
    return ISMIP7Atmosphere(esm=esm, scenario=scenario).available_years(var)


def clim_pool_years(esm, var):
    r"""Years the runs actually pool to build the reference climatology:
    historical + CLIM_SCENARIO, restricted to the CLIM_START-CLIM_END window.

    The window filter matters - both the control's ``compute_climatology`` and
    the projections' ``smb_scheme`` keep only years inside it, so an ESM with
    plenty of files outside the window still yields an empty pool. Deduplicated
    because a year present in both scenarios is one year of coverage.
    """
    return sorted({y for y in (atm_years(esm, "historical", var)
                               + atm_years(esm, CLIM_SCENARIO, var))
                   if CLIM_START <= y <= CLIM_END})


def clim_pool_gaps(esm, var, years=None):
    r"""Years of the climatology window the pool is MISSING.

    A partial pool is a failure, not a warning: the pooled mean is the baseline
    the aSMB anomalies are re-referenced to, so a core built from 15 of the 30
    years is referenced to a different mean than a sibling core with full
    coverage, while both are differenced against the same CTRL. That is the
    baseline mismatch this gate exists to catch, so it is held to the same
    standard as the run-scenario coverage check below.
    """
    have = clim_pool_years(esm, var) if years is None else years
    return sorted(set(range(CLIM_START, CLIM_END + 1)) - set(have))


def pool_status(esm, var, what_empty, what_partial):
    r"""``(bucket, detail)`` for the climatology pool: bucket is ``"ok"``,
    ``"partial"`` or ``"empty"``.

    The three are graded because the RUNTIME grades them: an empty pool makes
    the run take a different code path (the control raises, a projection
    silently falls back to full acabf(t)), while a partial pool runs and warns.
    A gate that blocked what the runtime happily runs would train people to
    ignore it, so a partial pool reports PARTIAL, not BLOCKED.

    The caller passes both consequences and the bucket picks between them, so
    the emptiness test lives in one place and the year list is read from disk
    once per (esm, var) rather than once per caller condition.
    """
    have = clim_pool_years(esm, var)
    gaps = clim_pool_gaps(esm, var, have)
    if not gaps:
        return "ok", None
    span = f"covers {have[0]}-{have[-1]}" if have else "no years"
    return ("partial" if have else "empty",
            f"{esm} {var} climatology pool (historical+{CLIM_SCENARIO}) "
            f"{span}, {len(gaps)} of {CLIM_START}-{CLIM_END} missing "
            f"({gaps[0]}..{gaps[-1]}): "
            f"{what_partial if have else what_empty}")


def ocx_atm_years():
    return ISMIP7Atmosphere(esm=OCX_ATMOSPHERE_SOURCE, scenario=OCX).available_years("acabf")


def ocx_ocean_cover():
    r"""The years both OCX ocean variables cover, or None."""
    ocean = ISMIP7Ocean(scenario=OCX, variant=_ocx_ocean())
    covers = [ocean.coverage(v) for v in ("tf", "so")]
    if None in covers:
        return None
    return max(c[0] for c in covers), min(c[1] for c in covers)


def ocean_cover(esm, scenario):
    r"""``(first, last)`` ocean forcing year on disk, as the reader sees it."""
    return ISMIP7Ocean(esm=esm, scenario=scenario).coverage("tf")


def msh_vertex_count(path):
    r"""The vertex count a gmsh ``.msh`` header gives, or None: the line after
    ``$Nodes`` holds it in version 2 files and as its second number in
    version 4."""
    try:
        with open(path, errors="replace") as f:
            for line in f:
                if line.startswith("$Nodes"):
                    parts = next(f).split()
                    return int(parts[1] if len(parts) >= 4 else parts[0])
    except (OSError, ValueError, StopIteration, IndexError):
        return None
    return None


def melt_calibration_missing(mesh_fn):
    r"""What stops a core from melting with its calibration.

    The resolver's own refusals come first: the file missing, a sidecar hash
    that does not match, a scaled K. Then the IMBIE2 basin grid the offsets
    are stamped through. Then the mesh: the tracked offsets were fitted on
    one mesh, the forward stamps them onto any mesh so a coarse probe runs,
    and this gate keeps a production core on the calibration's mesh. Naming
    a file with ISMIP7_DELTAT_PER_BASIN_NPZ (a refit on this mesh, or the
    tracked file itself) clears the mesh check. A legacy per-basin K named
    with ISMIP7_K_PER_BASIN_NPZ is checked for existence only."""
    import numpy as np
    try:
        npz = deltat_per_basin_npz()
        if npz is None:
            k_per_basin_npz()
            return []
    except (FileNotFoundError, ValueError) as e:
        return [f"melt calibration: {e}"]
    miss = []
    with np.load(npz) as data:
        recorded = str(data["imbie2_nc"]) if "imbie2_nc" in data else None
        K = float(data["K"])
    basins = imbie2_basin_path(recorded)
    if not os.path.exists(basins):
        miss.append(f"IMBIE2 basin grid ({basins}), which the melt "
                    f"calibration's offsets are stamped through")
    contract = melt_calibration_contract(npz)
    stem = os.path.splitext(os.path.basename(mesh_fn))[0]
    remedy = (f"Refit the offsets on this mesh (calibrate_deltaT.py --K {K:.3e}) "
              f"and name the file with ISMIP7_DELTAT_PER_BASIN_NPZ, or name the "
              f"tracked file there to run with its offsets")
    if contract and not os.environ.get("ISMIP7_DELTAT_PER_BASIN_NPZ"):
        if contract.get("mesh") != stem:
            miss.append(
                f"melt calibration: {os.path.basename(npz)} was fitted on "
                f"{contract.get('mesh')} and this core runs on {stem}. {remedy}")
        else:
            # Two builds of one mesh name differ cell by cell (IU's and
            # Rice's of the production mesh), and only the vertex count in
            # the header tells them apart.
            want = contract.get("vertices")
            have = msh_vertex_count(mesh_fn) if want else None
            if want and have is not None and have != int(want):
                miss.append(
                    f"melt calibration: {os.path.basename(npz)} was fitted on "
                    f"{contract.get('mesh_build', 'a build')} of {stem} "
                    f"({int(want)} vertices), and {mesh_fn} has {have}: "
                    f"another build of the same mesh. {remedy}")
    return miss


def shared_missing(warn=None):
    r"""Missing shared inputs. Non-fatal caveats are appended to ``warn``."""
    miss = []
    warn = warn if warn is not None else []
    # Front configuration: a mistyped law would otherwise surface only after
    # the forward's MAP load and initial solve.
    try:
        _calving_law()
        _calving_sigma_max()
    except ValueError as e:
        miss.append(str(e))
    # ISMIP7_MESH=checkpoint means the mesh inside the MAP; there is no file
    # to look for and the MAP check below covers it.
    mesh_fn = mesh_override() or mesh_filename(lc_coarse, lc, get_buffer_m())
    if not os.path.exists(mesh_fn):
        miss.append(f"mesh ({os.path.basename(mesh_fn)})")
    # The MAP the forward will actually load: ISMIP7_INVERSION if it names one
    # explicitly, else the one tagged with this geometry space, else the legacy
    # untagged (CG1) MAP it falls back to with a warning. A legacy MAP runs, but
    # its controls carry the CG1 front bias, so the run is a smoke test rather
    # than a result. An explicit override deliberately bypasses that lookup, so
    # setup_model raises on a missing file rather than falling back: report it
    # as a hard miss, exactly as the forward would.
    inv_override = os.environ.get("ISMIP7_INVERSION")
    inv = inv_override or os.path.join(MESH_DIR, map_basename(friction, lc))
    legacy = os.path.join(MESH_DIR, map_basename(friction, lc, geometry=False))
    if not os.path.exists(inv):
        if not inv_override and os.path.exists(legacy):
            warn.append(
                f"no {os.path.basename(inv)}; the forward would fall back to "
                f"{os.path.basename(legacy)} (inverted under a different "
                f"geometry space — smoke test only)"
            )
        else:
            miss.append(f"MAP ({os.path.basename(inv)})")
    # Same sidecar-name rule the solvers use: per-mesh preferred, shared file
    # as the fallback.
    bnd = sidecar_path(MESH_DIR, mesh_hint=mesh_fn)
    if not os.path.exists(bnd):
        miss.append(
            f"{os.path.basename(bnd)} (untracked — restore or regenerate "
            f"with make_boundary_ids.py)"
        )
    miss += melt_calibration_missing(mesh_fn)
    for d, pat, what in [
        (os.path.join(DATA_DIR, "bedmachine"), "*.nc", "BedMachine"),
        (os.path.join(DATA_DIR, "velocity"), "*.nc", "MEaSUREs velocity"),
    ]:
        if not glob.glob(os.path.join(d, pat)):
            miss.append(what)
    return miss


def racmo_ok():
    return os.path.exists(os.path.join(
        DATA_DIR, "racmo",
        "smbgl_monthlyS_ANT11_RACMO2.4p1_ERA5_197901_202312.nc",
    ))


def oi_ok():
    if root is None:
        return False
    return all(os.path.exists(_oi_climatology_path(root, v, oi_version))
               for v in ("tf", "so"))


def main():
    geom = _geometry_space()
    print(f"Preflight: lc={lc}, friction={friction}, geometry={geom}, "
          f"OI={oi_version}, climatology=historical+{CLIM_SCENARIO} "
          f"{CLIM_START}-{CLIM_END}")
    warn = []
    base_missing = shared_missing(warn)
    if base_missing:
        print(f"  SHARED inputs missing: {', '.join(base_missing)}")
    else:
        print("  Shared inputs (mesh, MAP, bndids, melt calibration, BedMachine, "
              "velocity): OK")
    for w in warn:
        print(f"  WARNING: {w}")
    print(f"  RACMO baseline: {'OK' if racmo_ok() else 'MISSING (acabf fallback)'}")
    print("  Status: READY = every input present; PARTIAL = the run proceeds "
          "but warns and its provenance is degraded; BLOCKED = missing input")
    print()

    for core, title, driver, esm, scenario in CORES:
        try:
            t_start, t_end = driver_period(driver, core)
        except ValueError as e:
            print(f"  core {core:2d}  BLOCKED  {title}  <- period: {e}")
            continue
        # the years the driver asks its readers for
        y0, y1 = int(math.floor(t_start + 1e-9)), forcing_year(t_end)
        miss = list(base_missing)
        degraded = []
        notes = []
        if scenario == "historical" or core == 11:
            # A cold start from the MAP. A start before 2015 backdates its
            # geometry by the Smith dH/dt (issue 117), which setup_model does
            # on the DG0 geometry only, and the drivers refuse a start before
            # the dH/dt window.
            try:
                if geometry_backdate_years(t_start) > 0.0 and geom != "dg0":
                    miss.append(f"a start at {t_start:g} backdates the geometry, "
                                f"which needs ISMIP7_GEOMETRY_SPACE=dg0")
            except ValueError as e:
                miss.append(f"geometry backdating: {e}")
        if core == 11 and _ocx_forcing() == "protocol":
            # projections/ocx.py refuses to start on anything less, so this
            # gate asks the same readers the same question.
            yrs = ocx_atm_years()
            gaps = sorted(set(range(y0, y1 + 1)) - set(yrs))
            if gaps:
                miss.append(f"OCX atmosphere acabf ({OCX_ATMOSPHERE_SOURCE}): "
                            + (f"{len(gaps)} of {y0}-{y1} missing" if yrs else "absent"))
            oc = ocx_ocean_cover()
            if oc is None or oc[0] > y0 or oc[1] < y1:
                miss.append(f"OCX ocean '{_ocx_ocean()}' tf/so"
                            + (f" covers {oc[0]}-{oc[1]}, need {y0}-{y1}" if oc else ": absent"))
            notes.append("K is fitted to the OI climatology, not to the OCX ocean: "
                         "read check_melt_bound.py --ocx first (discussion #48)")
        elif core == 11:
            notes.append("ISMIP7_OCX_FORCING=stopgap: RACMO2.4p1 + OI climatology, "
                         "not the ISMIP7 OCX product")
            if not racmo_ok():
                miss.append("RACMO (OCX SMB)")
            if not oi_ok():
                miss.append(f"OI climatology ({oi_version})")
        elif scenario is None:  # CTRL
            # the ESM's own ctrl ocean, which control/run.py refuses to start
            # without (icepack/ismip7#107)
            oc = ocean_cover(esm, "ctrl")
            if oc is None:
                miss.append(f"{esm}/ctrl ocean tf/so")
            elif oc[0] > y0 or oc[1] < y1 - 1:
                miss.append(f"ctrl ocean covers {oc[0]}-{oc[1]}, need {y0}-{y1}")
            if not racmo_ok():
                bucket, detail = pool_status(
                    esm, "acabf",
                    "the control would refuse to run",
                    "the constant SMB climatology would be a mean over part "
                    "of the window; the run warns and proceeds")
                if bucket == "empty":
                    miss.append(detail)
                elif bucket == "partial":
                    degraded.append(detail)
        else:
            # experiment.py re-references aSMB against the historical +
            # CLIM_SCENARIO acabf-anomaly pool and, on FileNotFoundError,
            # silently degrades to the full acabf(t) field - a different SMB
            # scheme from the one the CTRL this core is differenced against
            # uses. Only meaningful when RACMO is present: without it every
            # core degrades together, which the shared RACMO line reports.
            if racmo_ok():
                bucket, detail = pool_status(
                    esm, "acabf-anomaly",
                    "the run would silently fall back to full acabf(t), a "
                    "different SMB scheme than the CTRL",
                    "the aSMB re-reference baseline would differ from a "
                    "full-window sibling's; the run warns and proceeds")
                if bucket == "empty":
                    miss.append(detail)
                elif bucket == "partial":
                    degraded.append(detail)
            yrs = atm_years(esm, scenario) or atm_years(esm, scenario, "acabf")
            gaps = sorted(set(range(y0, y1 + 1)) - set(yrs))
            # The last year of a series may be absent: the reader persists
            # the last year on disk one year past the end (a CESM2-WACCM tree
            # fetched while the empty 2300 files were withdrawn stops at 2299,
            # discussion #8), so that is a note, not a missing input. Any
            # other gap is an error the reader raises on, so it blocks.
            bridged = gaps == [y1] and yrs and yrs[-1] == y1 - 1
            if not yrs:
                miss.append(f"{esm}/{scenario} atmosphere")
            elif bridged:
                notes.append(
                    f"atmosphere covers {yrs[0]}-{yrs[-1]}; {y1} is absent "
                    f"and the reader persists {y1 - 1} for it"
                )
            elif gaps:
                miss.append(
                    f"atmosphere covers {yrs[0]}-{yrs[-1]} with "
                    f"{len(gaps)} of {y0}-{y1} missing "
                    f"({gaps[0]}..{gaps[-1]})"
                )
            oc = ocean_cover(esm, scenario)
            if oc is None:
                miss.append(f"{esm}/{scenario} ocean tf/so")
            elif oc[0] > y0 or oc[1] < y1 - 1:
                miss.append(f"ocean covers {oc[0]}-{oc[1]}, need {y0}-{y1}")

        status = "BLOCKED" if miss else "PARTIAL" if degraded else "READY  "
        shown = miss + degraded + notes
        detail = "" if not shown else "  <- " + "; ".join(shown)
        print(f"  core {core:2d}  {status}  {title} {y0}-{y1}{detail}")


if __name__ == "__main__":
    main()
