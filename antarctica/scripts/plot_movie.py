#!/usr/bin/env python3
r"""Render the yearly checkpoints of a forward run into movie frames + mp4.

For each ``<prefix>_t<year>.h5`` checkpoint written by ``run_simulation`` (see
``ISMIP7_CHECKPOINT_EVERY_YR``), draw three panels on the run's own mesh:
thickness change since t=0 (``thickness - H_init``, the t=0 anchor every
checkpoint carries), surface speed, and thickness, with the t=0 ice extent
drawn as a black outline and the current extent in red so a moving calving
front is visible. A fourth panel plots the volume above flotation from the
run's ``*_timeseries.csv`` over the whole run, with the frame's own year
marked, so the maps can be read against the sea-level signal they produce.
The header carries the year plus the mass and VAF so a frame can be read on
its own.

Usage:
    python antarctica/scripts/plot_movie.py PREFIX [--out DIR] [--dh-max M]
        [--fps N] [--stride K]

PREFIX is the checkpoint basename without the ``_t<year>.h5`` suffix, e.g.
``antarctica/results/ctrl2015_cesm2_waccm_movie_velonly_32000``. Frames go to
``--out`` (default ``antarctica/figs/movie_<basename>/``) and the mp4 next
to them. Needs ``ffmpeg`` on PATH for the video; frames are still written
without it. Read-only with respect to the checkpoints; serial (one rank).
"""

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import firedrake as fd


def _year_of(path):
    m = re.search(r"_t(\d+(?:\.\d+)?)\.h5$", path)
    return float(m.group(1)) if m else None


def _load_timeseries(prefix):
    r"""``{year: (mass_gt, vaf_mm_sle)}`` and the same as sorted arrays."""
    fn = f"{prefix}_timeseries.csv"
    rows = {}
    if os.path.exists(fn):
        with open(fn) as f:
            for r in csv.DictReader(f):
                try:
                    rows[round(float(r["year"]), 1)] = (
                        float(r["mass_gt"]), float(r["vaf_mm_sle"]))
                except (KeyError, ValueError):
                    continue
    yrs = np.array(sorted(rows))
    return rows, (yrs, np.array([rows[y][1] for y in yrs]))


def _nearest(rows, year):
    if not rows:
        return None
    k = min(rows, key=lambda y: abs(y - year))
    return rows[k] if abs(k - year) <= 0.11 else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("prefix")
    ap.add_argument("--out", default=None)
    ap.add_argument("--dh-max", type=float, default=200.0,
                    help="symmetric colour limit for thickness change [m]")
    ap.add_argument("--fps", type=int, default=4)
    ap.add_argument("--stride", type=int, default=1,
                    help="use every K-th checkpoint")
    args = ap.parse_args()

    files = [f for f in glob.glob(f"{args.prefix}_t*.h5")
             if _year_of(f) is not None]
    files = sorted(files, key=_year_of)[:: args.stride]
    if not files:
        sys.exit(f"no checkpoints match {args.prefix}_t*.h5")

    base = os.path.basename(args.prefix)
    out = args.out or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "figs", f"movie_{base}")
    os.makedirs(out, exist_ok=True)
    ts, (ts_yr, ts_vaf) = _load_timeseries(args.prefix)

    mesh = None
    frame_paths = []
    outline0 = None
    for i, fn in enumerate(files):
        year = _year_of(fn)
        with fd.CheckpointFile(fn, "r") as chk:
            if mesh is None:
                mesh = chk.load_mesh()
            h = chk.load_function(mesh, name="thickness")
            h0 = chk.load_function(mesh, name="H_init")
            u = chk.load_function(mesh, name="velocity")
        Q = h.function_space()
        # Extent outlines: lumped CG1 projection of the cell indicator,
        # contoured at 0.5 (exact facets would need DMPlex; this is within
        # half a cell).
        Q1 = fd.FunctionSpace(mesh, "CG", 1)
        lump = fd.assemble(fd.TestFunction(Q1) * fd.dx).dat.data_ro

        def outline(field):
            ind = fd.Function(Q)
            ind.dat.data[:] = np.where(field.dat.data_ro > 1.0, 1.0, 0.0)
            f1 = fd.Function(Q1)
            f1.dat.data[:] = (fd.assemble(ind * fd.TestFunction(Q1) * fd.dx)
                              .dat.data_ro / lump)
            return f1

        if outline0 is None:
            outline0 = outline(h0)
        outline_now = outline(h)
        dh = fd.Function(Q, name="dh")
        dh.dat.data[:] = h.dat.data_ro - h0.dat.data_ro
        # Blank cells that never held ice so the buffer stays white.
        ice = (h.dat.data_ro > 0.5) | (h0.dat.data_ro > 0.5)
        dh.dat.data[~ice] = np.nan
        hm = fd.Function(Q)
        hm.dat.data[:] = np.where(ice, h.dat.data_ro, np.nan)
        speed = fd.Function(fd.FunctionSpace(mesh, "CG", 1))
        speed.dat.data[:] = np.maximum(
            np.hypot(u.dat.data_ro[:, 0], u.dat.data_ro[:, 1]), 1e-1)
        # Blank the buffer: the velocity there is the pinned ice-free
        # floor, not ice motion.
        speed.dat.data[outline_now.dat.data_ro < 0.5] = np.nan

        fig, axes = plt.subplots(1, 4, figsize=(26, 7),
                                 gridspec_kw={"width_ratios": [1, 1, 1, 0.75]})
        m = fd.tripcolor(dh, axes=axes[0], cmap="RdBu",
                         vmin=-args.dh_max, vmax=args.dh_max)
        fig.colorbar(m, ax=axes[0], shrink=0.7, label="h(t) - h(2015) [m]")
        axes[0].set_title("thickness change since t=0")
        m = fd.tripcolor(speed, axes=axes[1], cmap="viridis",
                         norm=LogNorm(vmin=1.0, vmax=5000.0))
        fig.colorbar(m, ax=axes[1], shrink=0.7, label="|u| [m/yr]")
        axes[1].set_title("surface speed")
        m = fd.tripcolor(hm, axes=axes[2], cmap="Blues", vmin=0, vmax=4000)
        fig.colorbar(m, ax=axes[2], shrink=0.7, label="h [m]")
        axes[2].set_title("thickness")
        # (4) Volume above flotation over the whole run, this frame marked.
        vax = axes[3]
        if len(ts_yr):
            dv = ts_vaf - ts_vaf[0]
            vax.plot(ts_yr, dv, color="0.6", lw=1.2)
            upto = ts_yr <= year + 1e-6
            vax.plot(ts_yr[upto], dv[upto], color="C0", lw=2.0)
            if upto.any():
                vax.plot(ts_yr[upto][-1], dv[upto][-1], "o", color="C0", ms=6)
            vax.set_xlim(ts_yr[0], ts_yr[-1])
            pad = max(0.05 * float(dv.max() - dv.min()), 1e-3)
            vax.set_ylim(float(dv.min()) - pad, float(dv.max()) + pad)
        vax.axhline(0.0, color="k", lw=0.6)
        vax.set_title("volume above flotation")
        vax.set_xlabel("year")
        vax.set_ylabel("VAF change since t=0 [mm SLE]")
        vax.grid(alpha=0.3)
        for ax in axes[:3]:
            fd.tricontour(outline0, axes=ax, levels=[0.5], colors="k",
                          linewidths=0.4)
            fd.tricontour(outline_now, axes=ax, levels=[0.5], colors="r",
                          linewidths=0.4)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
        head = f"{base}   t = {year:.1f}"
        mv = _nearest(ts, year)
        if mv is not None:
            head += f"   mass = {mv[0]/1e6:.4f} x 10^6 Gt   VAF = {mv[1]:.1f} mm SLE"
        fig.suptitle(head, fontsize=14)
        fig.tight_layout()
        fp = os.path.join(out, f"frame_{i:04d}.png")
        fig.savefig(fp, dpi=90)
        plt.close(fig)
        frame_paths.append(fp)
        print(f"  {os.path.basename(fn)} -> {os.path.basename(fp)}", flush=True)

    if shutil.which("ffmpeg") is None:
        print(f"{len(frame_paths)} frames in {out}; ffmpeg not found, no video")
        return
    mp4 = os.path.join(out, f"{base}.mp4")
    # -frames:v bounds the encode to the frames THIS run wrote: the image2
    # sequence otherwise reads frame_0000 upward to the first gap, so a
    # re-render with a larger stride would splice the tail of the old one.
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(args.fps),
           "-i", os.path.join(out, "frame_%04d.png"),
           "-frames:v", str(len(frame_paths)),
           "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", mp4]
    subprocess.run(cmd, check=True)
    print(f"{len(frame_paths)} frames -> {mp4}")


if __name__ == "__main__":
    main()
