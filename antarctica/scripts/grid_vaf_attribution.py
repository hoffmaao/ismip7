#!/usr/bin/env python
r"""Where the organisers' volume above flotation on the 8 km grid parts from
the model's, and how far shelf loss and grounding-line retreat move it
(issue #99). Read-only and serial.

    python antarctica/scripts/grid_vaf_attribution.py --out DIR \
        [--list FILE] [--no-perturb] [LABEL=CHECKPOINT ...]

A CHECKPOINT is a MAP or timing cache (``thickness``, ``bed``) or one year
of ISMIP7 annual output (``lithk``, ``topg``). ``--list`` names a file of
``LABEL CHECKPOINT`` lines. For each state this builds the writer's own
conservative overlap operator (``write_ismip7_output.overlap_operator``,
cached in DIR under a hash of the cell centroids, so the years of one run
share it), forms the pixel means the writer submits (``lithk`` over the whole
pixel, ``topg`` over the covered part) and sets ismip7-scalars' integrand
max(lithk - max(-topg, 0) rho_w/rho_i, 0) on them against the model's
max(h - h_f, 0) on grounded ice, taken to the grid by the same operator.
Every pixel is weighed by its map-plane area, with af2 = 1 and without
maxmask1, as in the issue 13 attribution and the grid VAF term of
``scalar_comparison_32km.md`` cause 2; the tool's own slvaf carries af2 and
maxmask1 on top of this (the area factor, cause 1 there, issue #97).
Pixels are classed as domain edge (partly covered by the mesh), grounding
line (covered and partly grounded) and interior. On a covered pixel grid
minus mesh splits exactly into

    convexity    = tool - max(mean lithk - mean h_f, 0)   >= 0
    cancellation = max(mean lithk - mean h_f, 0) - mesh   <= 0

the first because h_f is convex in the bed, the second because a pixel mean
lets the flotation deficit of its floating or ice-free part cancel the VAF of
its grounded part.

Then each state is perturbed. Floating ice thinned by 25, 50 and 75 % and
removed leaves the mesh's VAF where it was, so the grid's change is the
tool's error alone. The retreats remove the shelves and every marine
grounded cell less than X m of ice above flotation, for X from 25 to 800 m,
which moves the grounding line inland onto open water by one rule at every
resolution. ``delta_slvaf_mm`` is the tool's slvaf minus the model's, the
error a run would carry from that change.

Writes ``DIR/grid_vaf_<label>.json`` and ``DIR/pixels_<label>.npz`` (grid
minus mesh per pixel in Gt, as is and with the shelves removed, float32 on
the 761 x 761 grid). ``antarctica/reports/grid_vaf_resolution.md`` has the
September 2026 results; ``batch_runners/grid_vaf_attribution.script`` runs
it as a job.
"""
import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

from icepack2_tools.regrid import ISMIP7_DX, ISMIP7_NX, ISMIP7_NY  # noqa: E402

# the densities and ocean area of ismip7-scalars (compare_scalars.py)
RHO_I, RHO_W, RHO_FW, OCEAN_AREA = 917.0, 1024.0, 1000.0, 3.625e14
RATIO = RHO_W / RHO_I
PIXEL_AREA = ISMIP7_DX * ISMIP7_DX
GT_PER_MM = RHO_FW * OCEAN_AREA * 1.0e-3 / 1.0e12    # 362.5 Gt of VAF per mm of slvaf
GT_PER_PIXEL_M = PIXEL_AREA * RHO_I / 1.0e12         # 1 m of ice over one pixel, in Gt
TOL = 1.0e-9
THINNING = (0.25, 0.5, 0.75, 1.0)
RETREAT_HAF = (25.0, 50.0, 100.0, 200.0, 400.0, 800.0)
NAMES = (("thickness", "bed"), ("lithk", "topg"))


def evaluate(W, h, b, cell_area, keep=False):
    r"""Grid minus mesh VAF for one DG0 thickness and bed, by pixel class.

    ``W`` is the (pixels x cells) matrix of overlap areas. Grounded ice is
    the forward's: thicker than 1 m (the writer's ``lithk`` mask) and above
    flotation, which on DG0 is its ``s - s_float > 0``."""
    ice = h > 1.0
    hf = np.maximum(-b, 0.0) * RATIO
    grounded = ice & (h - hf > 0.0)
    lithk = np.where(ice, h, 0.0)
    mesh_cell = np.maximum(lithk - hf, 0.0) * grounded
    covered = W @ np.ones_like(h)
    cov = covered / PIXEL_AREA
    lithk_pix = (W @ lithk) / PIXEL_AREA
    topg_pix = np.zeros_like(cov)
    ok = covered > 0.0
    topg_pix[ok] = (W @ b)[ok] / covered[ok]
    g_pix = (W @ grounded.astype(float)) / PIXEL_AREA
    mesh_pix = (W @ mesh_cell) / PIXEL_AREA
    grid_pix = np.maximum(lithk_pix - np.maximum(-topg_pix, 0.0) * RATIO, 0.0)
    mean_haf = np.maximum(lithk_pix - (W @ hf) / PIXEL_AREA, 0.0)

    edge = (cov > TOL) & (cov < 1.0 - TOL)
    full = cov >= 1.0 - TOL
    gl = full & (g_pix > TOL) & (g_pix < 1.0 - TOL)
    interior = ~edge & ~gl
    d = (grid_pix - mesh_pix) * GT_PER_PIXEL_M
    convex = (grid_pix - mean_haf) * GT_PER_PIXEL_M * full
    cancel = (mean_haf - mesh_pix) * GT_PER_PIXEL_M * full
    mesh_gt = float((mesh_cell * cell_area).sum() * RHO_I / 1e12)
    out = {
        "mesh_vaf_gt": mesh_gt,
        "grid_vaf_gt": float(grid_pix.sum() * GT_PER_PIXEL_M),
        "grid_minus_mesh_gt": float(d.sum()),
        "edge_gt": float(d[edge].sum()),
        "gl_gt": float(d[gl].sum()),
        "interior_gt": float(d[interior].sum()),
        "gl_convexity_gt": float(convex[gl].sum()),
        "gl_cancellation_gt": float(cancel[gl].sum()),
        "interior_convexity_gt": float(convex[interior].sum()),
        "interior_cancellation_gt": float(cancel[interior].sum()),
        "gl_supply_gt": float(mesh_pix[gl].sum() * GT_PER_PIXEL_M),
        "gl_pixels": int(gl.sum()),
        "edge_pixels": int(edge.sum()),
        "edge_pixels_with_ice": int((edge & (lithk_pix > 0)).sum()),
        "grounded_area_km2": float(cell_area[grounded].sum() / 1e6),
        "floating_area_km2": float(cell_area[ice & ~grounded].sum() / 1e6),
        "floating_volume_km3": float((lithk * (ice & ~grounded) * cell_area).sum() / 1e9),
        # the regridded mesh integrand sums to the mesh's VAF when the grid
        # covers every cell, which the area gap in main() checks
        "mesh_sum_check": float(mesh_pix.sum() * GT_PER_PIXEL_M / max(mesh_gt, 1e-30) - 1.0),
    }
    if keep:
        out["_pixels"] = {"d_gt": d, "g_pix": g_pix, "cov": cov, "gl": gl}
        out["_cells"] = {"grounded": grounded, "ice": ice}
    return out


def perturbations(h, b, grounded, ice):
    r"""``(name, thickness, removed grounded cells)`` of every perturbed state."""
    floating = ice & ~grounded
    none = np.zeros_like(grounded)
    for f in THINNING:
        h2 = h.copy()
        h2[floating] *= 1.0 - f
        yield ("shelf_removed" if f == 1.0 else f"shelf_thinned_{int(f * 100)}"), h2, none
    haf = h - np.maximum(-b, 0.0) * RATIO
    for x in RETREAT_HAF:
        gone = grounded & (b < 0.0) & (haf < x)
        h2 = h.copy()
        h2[floating | gone] = 0.0
        yield f"retreat_haf{int(x)}", h2, gone


def cell_sizes(W, pixels, cell_area, ice):
    r"""Equilateral edge length, km, of the ice cells under ``pixels``."""
    sel = ((W.T @ pixels.astype(float)) > 0.0) & ice
    if not sel.any():
        return None
    L = np.sqrt(4.0 * cell_area[sel] / np.sqrt(3.0)) / 1e3
    return {"cells": int(sel.sum()), "median_km": float(np.median(L)),
            "area_weighted_mean_km": float((L * cell_area[sel]).sum() / cell_area[sel].sum())}


def _firedrake():
    r"""Firedrake, imported the way the writer imports it: icepack2 first."""
    import icepack2_tools.dual_friction  # noqa: F401  (icepack2 -> irksome import order)
    import firedrake
    return firedrake


def load_state(path):
    r"""The mesh, DG0 thickness and bed, and an annual file's ``sftgrf``."""
    fd = _firedrake()
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        Q = fd.FunctionSpace(mesh, "DG", 0)
        for hname, bname in NAMES:
            try:
                fh = chk.load_function(mesh, name=hname)
                fb = chk.load_function(mesh, name=bname)
            except Exception:  # noqa: BLE001  (not this file's names)
                continue
            break
        else:
            raise KeyError(f"{path}: neither thickness and bed nor lithk and topg")
        # the writer's topg is Function(Q).interpolate(b); on DG0 a copy
        h = fd.Function(Q).interpolate(fh).dat.data_ro.copy()
        b = fd.Function(Q).interpolate(fb).dat.data_ro.copy()
        sftgrf = (chk.load_function(mesh, name="sftgrf").dat.data_ro.copy()
                  if hname == "lithk" else None)
    return mesh, Q, h, b, sftgrf, [hname, bname]


def mesh_key(mesh):
    r"""A hash of the DG0 centroids in load order, which the operator's
    columns follow: two files of one mesh share it only if they load alike."""
    fd = _firedrake()
    X = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    return hashlib.sha1(np.ascontiguousarray(X).tobytes()).hexdigest()[:16]


def run_state(label, path, out, perturb=True):
    fd = _firedrake()
    from write_ismip7_output import overlap_operator
    t0 = time.time()
    mesh, Q, h, b, sftgrf, fields = load_state(path)
    key = mesh_key(mesh)
    cache = os.path.join(out, f"overlap_{len(h)}_{key}.npz")
    built = not os.path.exists(cache)
    t1 = time.time()
    W = overlap_operator(mesh, cache)
    t2 = time.time()
    cell_area = np.asarray(W.sum(axis=0)).ravel()
    true_area = fd.assemble(fd.TestFunction(Q) * fd.dx).dat.data_ro
    area_gap = float(np.max(np.abs(cell_area - true_area) / true_area))
    rec = {"label": label, "path": path, "cells": len(h), "fields": fields, "overlap": cache,
           "area_gap": area_gap, "load_s": t1 - t0, "overlap_s": t2 - t1, "experiments": {}}
    base = evaluate(W, h, b, cell_area, keep=True)
    px, cl = base.pop("_pixels"), base.pop("_cells")
    if sftgrf is not None:
        rec["sftgrf_mismatch_cells"] = int(((sftgrf > 0.5) != cl["grounded"]).sum())
    rec["cell_size_gl_pixels"] = cell_sizes(W, px["gl"], cell_area, cl["ice"])
    rec["experiments"]["as_is"] = base
    pix = {"as_is_d_gt": px["d_gt"], "as_is_g_pix": px["g_pix"], "cov": px["cov"]}
    if perturb:
        for name, h2, gone in perturbations(h, b, cl["grounded"], cl["ice"]):
            e = evaluate(W, h2, b, cell_area, keep=(name == "shelf_removed"))
            if name == "shelf_removed":
                pix["removed_d_gt"] = e.pop("_pixels")["d_gt"]
                e.pop("_cells")
            e["removed_grounded_area_km2"] = float(cell_area[gone].sum() / 1e6)
            e["delta_grid_minus_mesh_gt"] = e["grid_minus_mesh_gt"] - base["grid_minus_mesh_gt"]
            e["delta_mesh_vaf_gt"] = e["mesh_vaf_gt"] - base["mesh_vaf_gt"]
            e["delta_slvaf_mm"] = -e["delta_grid_minus_mesh_gt"] / GT_PER_MM
            rec["experiments"][name] = e
    np.savez_compressed(os.path.join(out, f"pixels_{label}.npz"),
                        **{k: v.reshape(ISMIP7_NY, ISMIP7_NX).astype("f4") for k, v in pix.items()})
    with open(os.path.join(out, f"grid_vaf_{label}.json"), "w") as f:
        json.dump(rec, f, indent=1)
    report(rec, built, time.time() - t0)
    return rec


def report(rec, built, seconds):
    cs = rec["cell_size_gl_pixels"] or {}
    print(f"[{rec['label']}] {rec['path']}\n    {rec['cells']} cells ({', '.join(rec['fields'])}), "
          f"loaded in {rec['load_s']:.0f} s, overlap {'built' if built else 'read'} in "
          f"{rec['overlap_s']:.0f} s, worst cell-area gap to the grid {rec['area_gap']:.1e}; "
          f"ice cells under grounding-line pixels: median {cs.get('median_km', float('nan')):.2f} km",
          flush=True)
    if "sftgrf_mismatch_cells" in rec:
        print(f"    cells where h > h_f disagrees with the file's sftgrf: {rec['sftgrf_mismatch_cells']}")
    print("    experiment          mesh VAF Gt  grid-mesh Gt    edge Gt      GL Gt  interior Gt"
          "  GL pixels  GL supply Gt  d(grid-mesh) Gt  slvaf mm  d(mesh) Gt  retreat km2")
    for name, e in rec["experiments"].items():
        print(f"    {name:18s} {e['mesh_vaf_gt']:12.1f} {e['grid_minus_mesh_gt']:+13.1f} "
              f"{e['edge_gt']:+10.1f} {e['gl_gt']:+10.1f} {e['interior_gt']:+12.1f} "
              f"{e['gl_pixels']:10d} {e['gl_supply_gt']:13.1f} "
              f"{e.get('delta_grid_minus_mesh_gt', 0.0):+16.1f} {e.get('delta_slvaf_mm', 0.0):+9.2f} "
              f"{e.get('delta_mesh_vaf_gt', 0.0):+11.1f} {e.get('removed_grounded_area_km2', 0.0):12.0f}")
    base = rec["experiments"]["as_is"]
    print(f"    mesh-sum identity {base['mesh_sum_check']:+.1e}; grounded {base['grounded_area_km2']:.4g} km2, "
          f"floating {base['floating_area_km2']:.4g} km2 and {base['floating_volume_km3']:.4g} km3; "
          f"{seconds:.0f} s", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, help="records, pixel planes and cached operators")
    p.add_argument("--list", help="a file of LABEL CHECKPOINT lines")
    p.add_argument("--no-perturb", action="store_true", help="the states as they are only")
    p.add_argument("states", nargs="*", metavar="LABEL=CHECKPOINT")
    a = p.parse_args()
    specs = []
    if a.list:
        with open(a.list) as f:
            specs += [tuple(line.split()) for line in f if line.strip() and not line.startswith("#")]
    specs += [tuple(s.split("=", 1)) for s in a.states]
    if not specs or any(len(s) != 2 for s in specs):
        p.error("name each state as LABEL=CHECKPOINT, or LABEL CHECKPOINT in --list")
    os.makedirs(a.out, exist_ok=True)
    for label, path in specs:
        run_state(label, path, a.out, perturb=not a.no_perturb)


if __name__ == "__main__":
    main()
