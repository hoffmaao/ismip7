#!/usr/bin/env python3
r"""Per-basin thermal-forcing offset deltaT_b at a fixed toolbox K.

The ISMIP7 ocean-forcing recommendation calibrates ONE dimensionless K for
the quadratic local melt law from the 4-term toolbox (K05/K50/K95) and then,
optionally, a per-basin offset deltaT_b of the thermal forcing so that each
IMBIE basin's present-day melt matches the remote-sensing total at that K
(`optimise_deltaT` in ismip7-antarctic-ocean-forcing/parameterisations/
parameter_selection_toolbox.py). This is the protocol's per-basin knob; a
per-basin K is not.

This script does that on OUR melt path: the same DG0 cell geometry, OI
climatology TF/so at the draft, constant mean-Antarctic slope and seawater
flotation test the forward melts with (calibrate_melt.forward_geometry), so
the deltaT the forward applies reproduces the basin totals by construction.

For each basin b and K:  M_b(dT) = sum_cells rho_i A_c m(TF_c + dT, S_c; K)
is monotone in dT (the law is TF |TF|), so the offset is a bracketed root of
M_b(dT) - M_obs_b on [-2, 2] K (the toolbox's search window); when no root
exists in the window the end point with the smaller residual is taken and
flagged, as the toolbox's grid argmin would.

Output: deltaT_per_basin_<lc>_K<K>.npz with basin_ids, deltaT_basin, K,
residual_gt (M_b(dT*) - M_obs_b), sensitivity_gt_per_K (dM_b/dT at dT*),
and the melt_slope / sin_alpha_ant / geometry_space provenance the forward
checks. Run with the same ISMIP7_* melt knobs as the forward.

    ISMIP7_LC=2000 python antarctica/scripts/calibrate_deltaT.py \
        [--K 8.5e-5 ...] [--out DIR]

The mesh comes from the MAP calibrate_melt.py names (ISMIP7_INV_H5 to name
another). Point a run at the result with ISMIP7_DELTAT_PER_BASIN_NPZ.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
RESULTS_DIR = os.path.join(os.path.dirname(HERE), "results")

import calibrate_melt as cm  # noqa: E402
from calibrate_melt import PETSc  # noqa: E402
from icepack2_tools.forcing import (  # noqa: E402
    quadratic_mixed_slope, sin_alpha_ant, melt_slope, _RHO_I, _K_PERCENTILES,
)
from icepack2_tools.runconfig import geometry_space  # noqa: E402

DT_WINDOW = (-2.0, 2.0)


def basin_totals(tf, sal, sin_a, K, floating, area, basin, bids, dT_b):
    r"""Gt/yr per basin at offsets ``dT_b`` (one per basin)."""
    out = np.zeros(len(bids))
    for i, bid in enumerate(bids):
        sel = floating & (basin == bid)
        if not sel.any():
            continue
        m = quadratic_mixed_slope(tf[sel] + dT_b[i], sal[sel], sin_a[sel], K=K)
        out[i] = float((m * area[sel]).sum()) * float(_RHO_I) / 1e12
    return out


def fit_deltaT(tf, sal, sin_a, K, floating, area, basin, bids, M_obs):
    from scipy.optimize import brentq
    dT = np.full(len(bids), np.nan)
    resid = np.full(len(bids), np.nan)
    sens = np.full(len(bids), np.nan)
    flagged = []
    for i, bid in enumerate(bids):
        sel = floating & (basin == bid)
        if not sel.any():
            flagged.append((bid, "no floating cells"))
            continue
        tf_b, s_b, a_b, A_b = tf[sel], sal[sel], sin_a[sel], area[sel]

        def M(d):
            m = quadratic_mixed_slope(tf_b + d, s_b, a_b, K=K)
            return float((m * A_b).sum()) * float(_RHO_I) / 1e12

        f = lambda d: M(d) - M_obs[i]  # noqa: E731
        lo, hi = DT_WINDOW
        if f(lo) * f(hi) < 0:
            dT[i] = brentq(f, lo, hi, xtol=1e-4)
        else:
            dT[i] = lo if abs(f(lo)) < abs(f(hi)) else hi
            flagged.append((bid, f"no root in [{lo:g}, {hi:g}] K, "
                                 f"took the end point"))
        resid[i] = f(dT[i])
        sens[i] = (M(dT[i] + 0.05) - M(dT[i] - 0.05)) / 0.1
    return dT, resid, sens, flagged


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--K", type=float, nargs="+", default=list(_K_PERCENTILES),
                    help="dimensionless K values (default: toolbox K05 K50 K95)")
    ap.add_argument("--out", default=RESULTS_DIR,
                    help="directory for deltaT_per_basin_<lc>_K<K>.npz")
    args = ap.parse_args()

    PETSc.Sys.Print("=== per-basin deltaT at fixed K (ISMIP7 toolbox optimise_deltaT, "
                    "on the forward's DG0 melt path) ===")
    cm._announce_obs_table()
    for p in (cm.INV_H5, cm.CLIM_TF, cm.CLIM_SO, cm.IMBIE2_NC, cm.OBS_CSV):
        if not os.path.exists(p):
            raise FileNotFoundError(p)
    bids, M_obs, sigma_obs = cm._load_obs()
    PETSc.Sys.Print(f"  Target: {float(M_obs.sum()):.1f} Gt/yr over {len(bids)} basins")

    mesh = cm._load_mesh()
    if cm.GEOMETRY != "dg0":
        raise SystemExit("run with ISMIP7_GEOMETRY_SPACE=dg0 (the forward's path)")
    g = cm.forward_geometry(mesh)
    mesh_x, mesh_y, draft = g["x"], g["y"], g["draft"]
    tf = cm._grid_interp(cm.CLIM_TF, "tf", mesh_x, mesh_y, draft=draft)
    sal = cm._grid_interp(cm.CLIM_SO, "so", mesh_x, mesh_y, draft=draft)
    basin = np.round(cm._grid_interp(cm.IMBIE2_NC, "basinNumber", mesh_x, mesh_y)).astype(int)
    sin_a = g["sin_a"]
    if cm.MELT_SLOPE != "ant":
        sin_a = np.minimum(sin_a, cm.SIN_ALPHA_CAP)
    floating, area = g["floating"], g["area"]
    PETSc.Sys.Print(f"  Slope: {cm.MELT_SLOPE}; floating cells {int(floating.sum())}; "
                    f"TF {tf[floating].min():.2f}..{tf[floating].max():.2f} K")

    os.makedirs(args.out, exist_ok=True)
    zero = np.zeros(len(bids))
    for K in args.K:
        M0 = basin_totals(tf, sal, sin_a, K, floating, area, basin, bids, zero)
        dT, resid, sens, flagged = fit_deltaT(tf, sal, sin_a, K, floating, area,
                                             basin, bids, M_obs)
        M1 = M_obs + resid
        PETSc.Sys.Print(f"\n  K = {K:.3e}: total {M0.sum():.0f} Gt/yr at dT=0, "
                        f"{M1.sum():.0f} with deltaT_b (obs {M_obs.sum():.0f})")
        PETSc.Sys.Print("  basin |  M_obs  sigma |  M(dT=0) |  dT* [K] |  M(dT*) | resid | dM/dT")
        for i, bid in enumerate(bids):
            PETSc.Sys.Print(f"  {bid:5d} | {M_obs[i]:6.1f} {sigma_obs[i]:5.1f} | "
                            f"{M0[i]:8.1f} | {dT[i]:+8.3f} | {M1[i]:7.1f} | "
                            f"{resid[i]:+5.1f} | {sens[i]:6.1f}")
        for bid, why in flagged:
            PETSc.Sys.Print(f"  [!] basin {bid}: {why}")
        fn = os.path.join(args.out, f"deltaT_per_basin_{cm.LC}_K{K:.3e}.npz")
        if mesh.comm.rank == 0:
            np.savez(fn, basin_ids=bids, deltaT_basin=dT, K=K,
                     M_obs=M_obs, M_dT0=M0, residual_gt=resid,
                     sensitivity_gt_per_K=sens,
                     melt_slope=melt_slope(), sin_alpha_ant=sin_alpha_ant(),
                     geometry_space=geometry_space(), obs_csv=cm.OBS_CSV,
                     imbie2_nc=cm.IMBIE2_NC, inversion=cm.INV_H5)
        PETSc.Sys.Print(f"  wrote {fn}")


if __name__ == "__main__":
    main()
