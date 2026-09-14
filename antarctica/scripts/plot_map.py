#!/usr/bin/env python
r"""Plot what an inversion derived: velocity, friction and fluidity.

    python antarctica/scripts/plot_map.py MAP.h5 [--out DIR] [--label NAME]
    python antarctica/scripts/plot_map.py A.h5 --diff B.h5 [--label A --label-b B]

Per MAP, six panels on the mesh: model speed, observed speed (MEaSUReS,
masked where no observation), their difference, the log friction adjustment
theta, the effective Weertman coefficient C = C_w0 exp(theta) on grounded
ice (C_w0 is the driving-stress anchor, rebuilt from the saved geometry and
observed velocity exactly as the forward does), and the log fluidity
adjustment phi (A = A_prior exp(phi)). With --diff, a second figure shows
B minus A for speed, theta and phi (same mesh required). A periodic
checkpoint carries no velocity: the speed panels are then left out.
"""
import argparse
import os
import re
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

from icepack2_tools.dual_friction import weertman_anchor            # noqa: E402  (icepack2 -> irksome first)
from icepack2_tools.grounding import height_above_flotation        # noqa: E402
from firedrake import CheckpointFile, Function, FunctionSpace       # noqa: E402

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
import matplotlib.tri as mtri                                       # noqa: E402
from matplotlib.colors import LogNorm, TwoSlopeNorm                 # noqa: E402


def load(map_path, m_slide=3.0):
    with CheckpointFile(map_path, "r") as chk:
        mesh = chk.load_mesh()
        f = {}
        for name in ("thickness", "surface", "bed", "log_friction", "log_fluidity",
                     "fluidity_prior", "velocity", "velocity_obs", "obs_mask"):
            try:
                f[name] = chk.load_function(mesh, name=name)
            except Exception:
                f[name] = None
        attrs = {k: chk.get_attr("/", k) for k in ("lc", "lc_coarse", "buffer_m", "mesh_basename")
                 if chk.has_attr("/", k)}
    xy = mesh.coordinates.dat.data_ro
    tri = mesh.coordinates.cell_node_map().values
    Q_g = f["thickness"].function_space()
    Q = f["log_friction"].function_space()
    H, s, b = f["thickness"], f["surface"], f["bed"]
    out = {"xy": xy / 1e3, "tri": tri, "attrs": attrs, "n_cells": tri.shape[0]}
    # cell-wise geometry (DG0 or CG1: interpolate to DG0 either way)
    Q0 = FunctionSpace(mesh, "DG", 0)
    out["haf"] = Function(Q0).interpolate(height_above_flotation(H, b)).dat.data_ro.copy()
    out["ice"] = Function(Q0).interpolate(H).dat.data_ro.copy() > 1.0
    theta = f["log_friction"].dat.data_ro.copy()
    phi = f["log_fluidity"].dat.data_ro.copy()
    out["theta"], out["phi"] = theta, phi
    out["nodal_controls"] = (theta.shape[0] == xy.shape[0])
    if f["velocity_obs"] is not None:
        C_w0 = weertman_anchor(H, s, f["velocity_obs"], m_slide, Q_g)
        c0 = Function(Q0).interpolate(C_w0).dat.data_ro.copy()
        th_cell = theta[tri].mean(axis=1) if out["nodal_controls"] else theta
        out["C"] = c0 * np.exp(th_cell)
        uo = f["velocity_obs"].dat.data_ro
        out["speed_obs"] = np.hypot(uo[:, 0], uo[:, 1])
        out["obs_mask"] = (f["obs_mask"].dat.data_ro.copy() > 0.5) if f["obs_mask"] is not None else np.ones(xy.shape[0], bool)
    if f["velocity"] is not None:
        u = f["velocity"].dat.data_ro
        out["speed"] = np.hypot(u[:, 0], u[:, 1])
    if f["fluidity_prior"] is not None:
        ap = f["fluidity_prior"].dat.data_ro.copy()
        out["A"] = ap * np.exp(phi) if ap.shape == phi.shape else None
    return out


def panel(ax, d, vals, title, cmap, norm=None, nodal=True, mask=None, units=""):
    tri = mtri.Triangulation(d["xy"][:, 0], d["xy"][:, 1], d["tri"])
    if mask is not None:
        tri.set_mask(~mask)
    if nodal:
        pc = ax.tripcolor(tri, vals, shading="gouraud", cmap=cmap, norm=norm, rasterized=True)
    else:
        pc = ax.tripcolor(tri, facecolors=vals, cmap=cmap, norm=norm, rasterized=True)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=10, loc="left")
    cb = plt.colorbar(pc, ax=ax, shrink=0.72, pad=0.02)
    cb.ax.tick_params(labelsize=7)
    if units:
        cb.set_label(units, fontsize=8)
    return pc


def figure_map(d, label):
    have_u = "speed" in d
    ncol = 3
    nrow = 2 if have_u else 1
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.3 * nrow + 0.4), dpi=110)
    axes = np.atleast_2d(axes)
    nodal = d["nodal_controls"]
    ice_nodes = None
    if have_u:
        sp = np.maximum(d["speed"], 0.1)
        so = np.maximum(d["speed_obs"], 0.1)
        panel(axes[0, 0], d, sp, "model speed |u|", "viridis", LogNorm(1, 4000), units="m/yr")
        panel(axes[0, 1], d, np.where(d["obs_mask"], so, np.nan), "observed speed (MEaSUReS)", "viridis", LogNorm(1, 4000), units="m/yr")
        diff = np.where(d["obs_mask"], d["speed"] - d["speed_obs"], np.nan)
        panel(axes[0, 2], d, diff, "model minus observed speed", "RdBu_r", TwoSlopeNorm(0, -200, 200), units="m/yr")
        r = 1
    else:
        r = 0
    panel(axes[r, 0], d, d["theta"], "log friction adjustment theta", "RdBu_r", TwoSlopeNorm(0, -3, 3), nodal=nodal)
    if "C" in d:
        grounded = d["ice"] & (d["haf"] > 0)
        C = np.where(grounded, np.maximum(d["C"], 1e-6), np.nan)
        panel(axes[r, 1], d, C, "C = C_w0 exp(theta), grounded", "plasma", LogNorm(1e-4, 1e-1), nodal=False, mask=grounded, units="MPa (m/yr)^-1/m")
    else:
        axes[r, 1].axis("off")
    panel(axes[r, 2], d, d["phi"], "log fluidity adjustment phi", "RdBu_r", TwoSlopeNorm(0, -3, 3), nodal=nodal)
    fig.suptitle(label, fontsize=12, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def figure_diff(a, b, la, lb):
    assert a["tri"].shape == b["tri"].shape, "different meshes"
    have_u = "speed" in a and "speed" in b
    ncol = 3 if have_u else 2
    fig, axes = plt.subplots(1, ncol, figsize=(4.6 * ncol, 4.7), dpi=110)
    k = 0
    if have_u:
        panel(axes[0], a, b["speed"] - a["speed"], f"speed: {lb} minus {la}", "RdBu_r", TwoSlopeNorm(0, -100, 100), units="m/yr"); k = 1
    panel(axes[k], a, b["theta"] - a["theta"], f"theta: {lb} minus {la}", "RdBu_r", TwoSlopeNorm(0, -2, 2), nodal=a["nodal_controls"])
    panel(axes[k + 1], a, b["phi"] - a["phi"], f"phi: {lb} minus {la}", "RdBu_r", TwoSlopeNorm(0, -2, 2), nodal=a["nodal_controls"])
    fig.suptitle(f"{lb} minus {la}", fontsize=12, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def slug(label):
    return re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")


def stats(d, label):
    ice = d["ice"]; grounded = ice & (d["haf"] > 0); floating = ice & (d["haf"] < 0)
    line = [f"{label}: {d['n_cells']} cells"]
    line.append(f"theta [{np.percentile(d['theta'], 1):+.2f}, {np.percentile(d['theta'], 99):+.2f}] (1-99%)")
    line.append(f"phi [{np.percentile(d['phi'], 1):+.2f}, {np.percentile(d['phi'], 99):+.2f}]")
    if "C" in d:
        Cg = d["C"][grounded]
        line.append(f"C grounded median {np.median(Cg):.3e}")
    if "speed" in d:
        m = d["obs_mask"]
        err = np.abs(d["speed"] - d["speed_obs"])[m]
        line.append(f"|u|-|u_obs| median {np.median(err):.1f} m/yr, p90 {np.percentile(err, 90):.1f}")
    return "  ".join(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("map"); ap.add_argument("--diff", default=None)
    ap.add_argument("--label", default=None); ap.add_argument("--label-b", default=None)
    ap.add_argument("--out", default=os.path.join(_ROOT, "..", "figs", "maps"))
    ap.add_argument("--m-slide", type=float, default=float(os.environ.get("ISMIP7_M_SLIDE", "3.0")))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    la = a.label or os.path.basename(a.map).replace("inversion_icepack2_", "").replace(".h5", "")
    d = load(a.map, a.m_slide)
    print(stats(d, la))
    fig = figure_map(d, la); p = os.path.join(a.out, f"map_{slug(la)}.png"); fig.savefig(p); plt.close(fig); print("wrote", p)
    if a.diff:
        lb = a.label_b or os.path.basename(a.diff).replace("inversion_icepack2_", "").replace(".h5", "")
        e = load(a.diff, a.m_slide); print(stats(e, lb))
        fig = figure_map(e, lb); p = os.path.join(a.out, f"map_{slug(lb)}.png"); fig.savefig(p); plt.close(fig); print("wrote", p)
        fig = figure_diff(d, e, la, lb); p = os.path.join(a.out, f"diff_{slug(lb)}_minus_{slug(la)}.png"); fig.savefig(p); plt.close(fig); print("wrote", p)


if __name__ == "__main__":
    main()
