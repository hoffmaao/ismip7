#!/usr/bin/env python3
"""Observed calving flux at BedMachine's own 500 m ice-shelf front.

Every edge between a floating pixel (mask 3) and an ocean pixel (mask 0)
carries H * max(u . n, 0) * 500 m, H the floating pixel's BedMachine
thickness, n the outward normal, u MEaSUREs 450 m v2 at the nearest pixel to
a point `inset` pixels inward of the floating pixel (inset 0 is the pixel
itself; insets avoid velocity pixels that mix ice and water). The reference
for a model front's flux (front_flux_check.py):

    python antarctica/scripts/bm_front_flux.py BEDMACHINE.nc MEASURES_450m_v2.nc

"""
import sys
import numpy as np
import netCDF4

bm = netCDF4.Dataset(sys.argv[1]); vel = netCDF4.Dataset(sys.argv[2])
x = np.asarray(bm["x"][:], float); y = np.asarray(bm["y"][:], float)
mask = np.asarray(bm["mask"][:]); H = np.asarray(bm["thickness"][:], float)
dx = abs(x[1] - x[0]); sx = np.sign(x[1] - x[0]); sy = np.sign(y[1] - y[0])
vxg = np.asarray(vel["x"][:], float); vyg = np.asarray(vel["y"][:], float)
VX = np.asarray(vel["VX"][:], float); VY = np.asarray(vel["VY"][:], float)
VX[~np.isfinite(VX) | (np.abs(VX) > 1e5)] = np.nan
VY[~np.isfinite(VY) | (np.abs(VY) > 1e5)] = np.nan
fl = mask == 3; oc = mask == 0
ny_, nx_ = mask.shape
edges = []
for di, dj in ((0, 1), (0, -1), (1, 0), (-1, 0)):
    ii, jj = np.nonzero(fl)
    ti, tj = ii + di, jj + dj
    inside = (ti >= 0) & (ti < ny_) & (tj >= 0) & (tj < nx_)
    ii, jj, ti, tj = ii[inside], jj[inside], ti[inside], tj[inside]
    hit = oc[ti, tj]
    edges.append((ii[hit], jj[hit], di, dj))
n_all = sum(len(e[0]) for e in edges)
print(f"BedMachine floating/ocean edges: {n_all} ({n_all * dx / 1e3:,.0f} km)")
hs = np.concatenate([H[e[0], e[1]] for e in edges])
print(f"front-pixel thickness mean {hs.mean():.1f} m, median {np.median(hs):.1f} m, p10 {np.percentile(hs,10):.1f}, p90 {np.percentile(hs,90):.1f}")
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
