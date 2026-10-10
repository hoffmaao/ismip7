#!/usr/bin/env python3
"""Observed calving flux at a 500 m ice front.

Every edge between a front pixel and an ocean pixel carries
H * max(u . n, 0) * 500 m, H the front pixel's BedMachine thickness, n the
outward normal, u MEaSUREs 450 m v2 at the nearest pixel to a point `inset`
pixels inward of the front pixel (inset 0 is the pixel itself; insets avoid
velocity pixels that mix ice and water). The reference for a model front's
flux (front_flux_check.py).

By default the front is BedMachine's own: floating pixels (mask 3) beside
ocean (mask 0). With --greene-tif it is the marine front of that Greene et
al. (2022) ice mask (obs_icemask.icemask_tif): ice pixels beside marine
ones (obs_icemask.classify), the front the _front<year> meshes follow
(issue #167); an ice pixel BedMachine holds no ice on is counted and carries
nothing.

    python antarctica/scripts/bm_front_flux.py BEDMACHINE.nc MEASURES_450m_v2.nc
    python antarctica/scripts/bm_front_flux.py BEDMACHINE.nc MEASURES_450m_v2.nc \\
        --greene-tif icemask_greene_2015_AntarcticaObsISMIP7-v1.2.tif

"""
import argparse
import os
import sys

import numpy as np
import netCDF4

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("bedmachine")
parser.add_argument("velocity")
parser.add_argument("--greene-tif", default=None)
args = parser.parse_args()

if args.greene_tif:
    import rasterio
    from icepack2_tools.obs_icemask import ICE, MARINE, _aligned_window, read_classes
    classes, tr = read_classes(args.greene_tif, args.bedmachine)
    with rasterio.open(f"netcdf:{args.bedmachine}:thickness") as src:
        bounds = rasterio.transform.array_bounds(*classes.shape, tr)
        H = src.read(1, window=_aligned_window(
            src, (bounds[0], bounds[1], bounds[2], bounds[3]))).astype(float)
    x = tr.c + tr.a * (np.arange(classes.shape[1]) + 0.5)
    y = tr.f + tr.e * (np.arange(classes.shape[0]) + 0.5)
    front, water = classes == ICE, classes == MARINE
    label = f"Greene {os.path.basename(args.greene_tif)} ice/marine"
else:
    bm = netCDF4.Dataset(args.bedmachine)
    x = np.asarray(bm["x"][:], float); y = np.asarray(bm["y"][:], float)
    mask = np.asarray(bm["mask"][:]); H = np.asarray(bm["thickness"][:], float)
    front, water = mask == 3, mask == 0
    label = "BedMachine floating/ocean"

vel = netCDF4.Dataset(args.velocity)
dx = abs(x[1] - x[0]); sx = np.sign(x[1] - x[0]); sy = np.sign(y[1] - y[0])
vxg = np.asarray(vel["x"][:], float); vyg = np.asarray(vel["y"][:], float)
VX = np.asarray(vel["VX"][:], float); VY = np.asarray(vel["VY"][:], float)
VX[~np.isfinite(VX) | (np.abs(VX) > 1e5)] = np.nan
VY[~np.isfinite(VY) | (np.abs(VY) > 1e5)] = np.nan
ny_, nx_ = front.shape
edges = []
for di, dj in ((0, 1), (0, -1), (1, 0), (-1, 0)):
    ii, jj = np.nonzero(front)
    ti, tj = ii + di, jj + dj
    inside = (ti >= 0) & (ti < ny_) & (tj >= 0) & (tj < nx_)
    ii, jj, ti, tj = ii[inside], jj[inside], ti[inside], tj[inside]
    hit = water[ti, tj]
    edges.append((ii[hit], jj[hit], di, dj))
n_all = sum(len(e[0]) for e in edges)
print(f"{label} edges: {n_all} ({n_all * dx / 1e3:,.0f} km)")
hs = np.concatenate([H[e[0], e[1]] for e in edges])
print(f"front-pixel thickness mean {hs.mean():.1f} m, median {np.median(hs):.1f} m, p10 {np.percentile(hs,10):.1f}, p90 {np.percentile(hs,90):.1f}")
if args.greene_tif:
    held = hs > 0
    print(f"front pixels BedMachine holds ice on: {held.mean():.1%}; their thickness "
          f"mean {hs[held].mean():.1f} m, median {np.median(hs[held]):.1f} m")
for inset in (0, 1, 2, 4, 8):
    tot = 0.0; n_ok = 0; un_sum = 0.0
    for ii, jj, di, dj in edges:
        si = np.clip(ii - di * inset, 0, ny_ - 1); sj = np.clip(jj - dj * inset, 0, nx_ - 1)
        px, py = x[sj], y[si]
        ci = np.clip(np.rint((px - vxg[0]) / (vxg[1] - vxg[0])).astype(int), 0, len(vxg) - 1)
        ri = np.clip(np.rint((py - vyg[0]) / (vyg[1] - vyg[0])).astype(int), 0, len(vyg) - 1)
        u = VX[ri, ci]; v = VY[ri, ci]
        un = u * dj * sx + v * di * sy
        ok = np.isfinite(un)
        tot += float((H[ii, jj][ok] * np.maximum(un[ok], 0.0)).sum() * dx)
        n_ok += int(ok.sum()); un_sum += float(un[ok].sum())
    gt = tot * 917.0 / 1e12
    print(f"inset {inset} px: velocity on {n_ok / n_all:.1%} of edges, mean u.n {un_sum / max(n_ok,1):.1f} m/yr, "
          f"flux {gt:,.1f} Gt/yr ({gt * n_all / max(n_ok, 1):,.1f} if the edges without velocity carried the same)")
