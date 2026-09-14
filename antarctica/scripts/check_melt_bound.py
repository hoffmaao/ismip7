#!/usr/bin/env python
r"""Where does the melt parameterisation exceed the variable request's bound?

    python antarctica/scripts/check_melt_bound.py [--npz results/calibrated_K_per_basin_2000.npz]

The ISMIP7 variable request gives ``libmassbffl`` an AIS minimum of
-0.008 kg m-2 s-1 with severity ``error``. In ice-equivalent thickness that is
275.3 m/yr, and a 10-year Ua-mesh ssp585 reached -0.0117 (402.6 m/yr) on
grounding-zone cells. Two readings fit that: the parameterisation is too strong
somewhere, or the writer's ``no_floating_ice`` fill policy reports one hot cell
as the whole 8 km pixel's value, which the request's own convention asks for.

This separates them on the model side, before any regridding. It rebuilds the
melt field the forward applies at the reference geometry, with the calibrated
per-basin K, and reports how much floating AREA sits past the bound. A handful
of small cells says the grid value is a fill-policy artefact of a nearly
ice-free pixel; a broad region says the melt is too strong.

Scope: the reference state, with the OI thermal-forcing climatology and the
BedMachine geometry. A projection's thermal forcing warms above the
climatology and its shelves thin, so a clean result here bounds the
parameterisation itself rather than what any particular run reaches.

Serial. Reuses calibrate_melt's loaders, so it needs the same inputs: a MAP for
the mesh, the OI climatology, the IMBIE2 basins and BedMachine.

Measured on the Ua 2 km mesh, 14 September 2026, with the per-basin K
calibrated against the re-released observation table: maximum 71.1 m/yr, 99th
percentile 22.2 m/yr, area mean 0.77 m/yr over 1 512 899 km2 of floating ice,
and nothing at all past the 275.3 m/yr bound.
"""
import argparse
import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))
sys.path.insert(0, _ROOT)

import numpy as np                                                    # noqa: E402
import firedrake as fd                                                # noqa: E402
from firedrake import FunctionSpace, assemble, dx                     # noqa: E402
from firedrake.petsc import PETSc                                     # noqa: E402

import calibrate_melt as cm                                           # noqa: E402
from icepack2_tools.forcing import quadratic_mixed_slope              # noqa: E402
# The same year and density the writer converts with, so the bound compared
# here is the one the checker applies.
from icepack2_tools.ismip7_output import RHO_I, SECONDS_PER_YEAR      # noqa: E402

# ISMIP7_variable_request.csv, min_value_ais for libmassbffl, severity error.
BOUND_KG_M2_S = -0.008


def m_per_yr(kg_m2_s):
    r"""kg m-2 s-1 of ice to m/yr of ice thickness."""
    return kg_m2_s * SECONDS_PER_YEAR / RHO_I


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", default=None,
                    help="calibrated_K_per_basin_<lc>.npz (default: the one "
                         "beside this run's lc)")
    ap.add_argument("--bound", type=float, default=BOUND_KG_M2_S,
                    help="the request's min_value_ais, kg m-2 s-1")
    a = ap.parse_args()

    bound_m_yr = abs(m_per_yr(a.bound))
    PETSc.Sys.Print(f"=== melt against the request bound {a.bound} kg m-2 s-1 "
                    f"({bound_m_yr:.1f} m/yr ice) ===")

    npz_path = a.npz or os.path.join(
        os.path.dirname(os.path.dirname(_ROOT)), "antarctica", "results",
        f"calibrated_K_per_basin_{cm.LC}.npz")
    d = np.load(npz_path, allow_pickle=True)
    PETSc.Sys.Print(f"  K from {os.path.basename(npz_path)}"
                    + (f", calibrated against {d['obs_csv']}"
                       if "obs_csv" in d else ""))

    mesh = cm._load_mesh()
    Q = FunctionSpace(mesh, "CG", 1)
    PETSc.Sys.Print(f"  Mesh: {mesh.num_vertices()} vertices, "
                    f"{mesh.num_cells()} cells")

    bed, thk, sur, msk = cm._interp_bedmachine(mesh, Q)
    h_np = thk.dat.data_ro
    s_np = sur.dat.data_ro
    mask_np = msk.dat.data_ro
    x = mesh.coordinates.dat.data_ro[:, 0]
    y = mesh.coordinates.dat.data_ro[:, 1]

    draft = np.minimum(s_np - h_np, 0.0)
    tf = cm._grid_interp(cm.CLIM_TF, "tf", x, y, draft=draft)
    sal = cm._grid_interp(cm.CLIM_SO, "so", x, y, draft=draft)
    sin_a = np.minimum(cm._compute_sin_alpha(mesh, thk, sur), cm.SIN_ALPHA_CAP)
    floating = (np.round(mask_np).astype(int) == 3)

    # The per-basin K field the forward stamps onto the mesh. K_field in the
    # npz was built on the calibration mesh; rebuild it here from K_basin so
    # this runs against any mesh.
    basin = np.round(cm._grid_interp(cm.IMBIE2_NC, "basinNumber", x, y)).astype(int)
    K = np.zeros_like(tf)
    for bid, kb in zip(d["basin_ids"], d["K_basin"]):
        if np.isfinite(kb):
            K[basin == int(bid)] = kb

    melt = np.where(floating, quadratic_mixed_slope(tf, sal, sin_a, K=K), 0.0)

    v = fd.TestFunction(Q)
    area = assemble(v * dx).dat.data_ro                 # nodal area weights, m^2
    afl = float(area[floating].sum())

    over = floating & (melt > bound_m_yr)
    a_over = float(area[over].sum())
    PETSc.Sys.Print(
        f"\n  floating area           {afl / 1e6:12.1f} km^2\n"
        f"  melt max                {melt.max():12.1f} m/yr\n"
        f"  melt p99 (floating)     {np.quantile(melt[floating], 0.99):12.1f} m/yr\n"
        f"  melt area-mean          "
        f"{float((melt * area)[floating].sum()) / afl:12.2f} m/yr\n"
        f"  nodes past the bound    {int(over.sum()):12d} of {int(floating.sum())}\n"
        f"  area past the bound     {a_over / 1e6:12.1f} km^2 "
        f"({100 * a_over / afl:.3f}% of floating)")

    if over.any():
        PETSc.Sys.Print("\n  worst nodes (x km, y km, melt m/yr, TF K, draft m, "
                        "sin_alpha, area km^2):")
        idx = np.argsort(-melt)[:10]
        for i in idx:
            PETSc.Sys.Print(
                f"    {x[i] / 1e3:9.1f} {y[i] / 1e3:9.1f} {melt[i]:9.1f} "
                f"{tf[i]:6.2f} {draft[i]:8.1f} {sin_a[i]:9.2e} "
                f"{area[i] / 1e6:8.2f}")
        # An 8 km pixel is 64 km^2. A node whose own area is a small fraction
        # of that cannot fill a pixel on its own, so its value reaching the
        # grid means the pixel carried little other floating ice.
        PETSc.Sys.Print(
            f"\n  median area of a node past the bound: "
            f"{np.median(area[over]) / 1e6:.2f} km^2, against 64 km^2 for an "
            f"8 km pixel")
    else:
        PETSc.Sys.Print("\n  nothing past the bound on this mesh")


if __name__ == "__main__":
    main()
