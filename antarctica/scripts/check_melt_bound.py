#!/usr/bin/env python
r"""How far does the draft-slope cap move melt past the variable request's bound?

    python antarctica/scripts/check_melt_bound.py [--npz results/calibrated_K_per_basin_2000.npz]

The ISMIP7 variable request gives ``libmassbffl`` an AIS minimum of
-0.008 kg m-2 s-1 with severity ``error``. In ice-equivalent thickness that is
275.3 m/yr, and a 10-year Ua-mesh ssp585 reached -0.0117 (402.6 m/yr) on
grounding-zone cells. Two readings fit that: the parameterisation is too strong
somewhere, or the writer's ``no_floating_ice`` fill policy reports one hot cell
as the whole 8 km pixel's value, which the request's own convention asks for.

This looks at the first on the model side, before any regridding. At the
reference geometry, with the calibrated per-basin K, it evaluates the melt four
times, capped and uncapped for each of two halves, and for each row reports
the maximum, the 99th percentile, the area mean, the integrated total against
the ``obs_total_gtyr`` the K was fitted to, and how much floating AREA sits
past the bound.

The calibration half reproduces calibrate_melt on CG1 nodes: BedMachine
interpolated with its raster ``surface`` and ``mask``, and grad(draft)
projected onto CG1.

* capped: sin_alpha capped at ISMIP7_SIN_ALPHA_CAP (calibrate_melt's default
  5e-3), the slope K was fitted against;
* uncapped: the same slope with no cap.

The forward half reproduces the forward on DG0 cells, the field the forward
melts with under DG0 geometry: bed and thickness sampled onto the cells, the
surface from flotation, s = max(b + H, (1 - 917/1024) H), as simulation.py
builds it, ``forcing.compute_sin_alpha``'s cell slope, thermal forcing and
salinity at each cell centroid and its own draft, and the forward callback's
``haf <= 0`` floating test.

* uncapped: as the forward runs today;
* capped: the cell slope capped at the same value, which is where a cap inside
  ``forcing.compute_sin_alpha`` would act.

The two halves use different floating masks and different quadrature, nodal
area weights against cell areas, so their totals compare in magnitude and
differ in detail.

The bound is ``min_value_ais`` for ``libmassbffl`` in the same bundled request
table the writer reads, converted with the writer's year and ice density, so it
is the bound the compliance checker applies.

Scope: the reference state, with the OI thermal-forcing climatology and the
BedMachine geometry. A projection's thermal forcing warms above the
climatology and its shelves thin.

Serial. Reuses calibrate_melt's loaders, so it needs the same inputs: a MAP for
the mesh, the OI climatology, the IMBIE2 basins and BedMachine.

Measured on the Ua 2 km mesh, September 2026, with
calibrated_K_per_basin_2000.npz calibrated against the re-released observation
table:

* calibration half, capped at 5e-3, over 1 512 899 km2 and 47 288 nodes:
  maximum 71.1 m/yr, 99th percentile 22.2, area mean 0.77, 1067 Gt/yr, which
  is the 1067.4 Gt/yr target K was fitted to, and nothing past the bound;
* calibration half, uncapped: maximum 1804.9 m/yr, 99th percentile 256.3, area
  mean 4.18, 5803 Gt/yr, and 421 nodes past the bound over 2920.6 km2 with a
  median node area of 6.33 km2;
* forward half, uncapped, over 1 631 466 km2 and 85 820 cells: maximum
  1144.1 m/yr, 99th percentile 59.0, area mean 1.16, 1732 Gt/yr, and 97 cells
  past the bound over 388.2 km2 with a median cell area of 3.61 km2;
* forward half, capped at 5e-3: maximum 57.8 m/yr, 99th percentile 13.0, area
  mean 0.43, 646 Gt/yr, and nothing past the bound.

The uncapped forward half integrates within 7% of the 1860 Gt/yr the 10-year
run booked, which validates it against a real run. As the two halves stand, the
forward melts 1.62 times what K was fitted to. Capping the forward's slope
undershoots the target by 39%, so neither convention on its own reconciles the
halves, which also differ in floating area, mask and quadrature. Recalibrating
K through the forward's own cell by cell melt path, under whichever slope
convention is chosen, reconciles them by construction.

An earlier form of this script lifted the forward's slope onto CG1 nodes and
melted it with CG1 forcing and the raster mask, so its forward rows,
4293 Gt/yr uncapped and 1028 Gt/yr capped, reproduced neither half and are
superseded.
"""
import argparse
import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))
sys.path.insert(0, _ROOT)

import numpy as np                                                    # noqa: E402
import rasterio                                                       # noqa: E402
import firedrake as fd                                                # noqa: E402
from firedrake import (FunctionSpace, VectorFunctionSpace,            # noqa: E402
                       assemble, dx)
from firedrake.petsc import PETSc                                     # noqa: E402

import calibrate_melt as cm                                           # noqa: E402
from icepack2_tools.forcing import (quadratic_mixed_slope,            # noqa: E402
                                    compute_sin_alpha,
                                    _RHO_ICE, _RHO_WATER)
from icepack2_tools.geometry import sample_to_geometry                # noqa: E402
from icepack2_tools.runconfig import raster_sample                    # noqa: E402
# The same year and density the writer converts with, so the bound compared
# here is the one the checker applies.
from icepack2_tools.ismip7_output import RHO_I, SECONDS_PER_YEAR      # noqa: E402
from icepack2_tools.regrid import ISMIP7_DX                           # noqa: E402
from write_ismip7_output import request_table                         # noqa: E402

# The ice to seawater density ratio simulation.py builds the surface with.
RHO_RATIO = 917.0 / 1024.0


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
    a = ap.parse_args()

    bound = float(request_table()["libmassbffl"]["min_value_ais"])
    bound_m_yr = abs(m_per_yr(bound))
    PETSc.Sys.Print(f"=== melt against the request bound {bound} kg m-2 s-1 "
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
    Q_g = FunctionSpace(mesh, "DG", 0)
    PETSc.Sys.Print(f"  Mesh: {mesh.num_vertices()} vertices, "
                    f"{mesh.num_cells()} cells")

    def k_at(xs, ys):
        r"""The per-basin K the forward stamps onto the mesh. K_field in the
        npz was built on the calibration mesh; rebuild it here from K_basin
        so this runs against any mesh."""
        basin = np.round(
            cm._grid_interp(cm.IMBIE2_NC, "basinNumber", xs, ys)).astype(int)
        K = np.zeros(len(xs))
        for bid, kb in zip(d["basin_ids"], d["K_basin"]):
            if np.isfinite(kb):
                K[basin == int(bid)] = kb
        return K

    def half(xs, ys, draft, sin_a, floating, area):
        r"""One half's melt inputs, all aligned with its own dofs."""
        return {"x": xs, "y": ys, "draft": draft, "sin_a": sin_a,
                "floating": floating, "area": area,
                "tf": cm._grid_interp(cm.CLIM_TF, "tf", xs, ys, draft=draft),
                "sal": cm._grid_interp(cm.CLIM_SO, "so", xs, ys, draft=draft),
                "K": k_at(xs, ys)}

    # The calibration half, as calibrate_melt builds it: BedMachine
    # interpolated onto CG1 nodes with its raster surface and mask, and
    # grad(draft) projected onto CG1.
    bed, thk, sur, msk = cm._interp_bedmachine(mesh, Q)
    x = mesh.coordinates.dat.data_ro[:, 0]
    y = mesh.coordinates.dat.data_ro[:, 1]
    calibration = half(
        x, y, np.minimum(sur.dat.data_ro - thk.dat.data_ro, 0.0),
        cm._compute_sin_alpha(mesh, thk, sur),
        np.round(msk.dat.data_ro).astype(int) == 3,
        assemble(fd.TestFunction(Q) * dx).dat.data_ro)

    # The forward half, as the forward melts cell by cell under DG0 geometry:
    # bed and thickness sampled onto the cells, the surface from flotation as
    # simulation.py builds it, forcing.compute_sin_alpha's DG0 slope, forcing
    # at each cell centroid and its own draft, and the callback's haf <= 0
    # floating test.
    bm = cm._bedmachine_path()
    b_dg = sample_to_geometry(rasterio.open(f"netcdf:{bm}:bed"), Q_g, Q,
                              method=raster_sample())
    h_dg = sample_to_geometry(rasterio.open(f"netcdf:{bm}:thickness"), Q_g, Q,
                              method=raster_sample())
    s_dg = fd.Function(Q_g).interpolate(
        fd.max_value(b_dg + h_dg, (1.0 - RHO_RATIO) * h_dg))
    xy_dg = fd.Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    b_np, h_np, s_np = b_dg.dat.data_ro, h_dg.dat.data_ro, s_dg.dat.data_ro
    haf = s_np - (b_np + (_RHO_WATER / _RHO_ICE) * np.maximum(-b_np, 0.0))
    forward = half(
        xy_dg[:, 0], xy_dg[:, 1], np.minimum(s_np - h_np, 0.0),
        compute_sin_alpha({"Q": Q, "V": VectorFunctionSpace(mesh, "CG", 1),
                           "Q_g": Q_g, "h": h_dg, "s": s_dg}),
        haf <= 0,
        assemble(fd.TestFunction(Q_g) * dx).dat.data_ro)

    cap = cm.SIN_ALPHA_CAP
    cases = [
        (f"calibration half, capped at {cap:.0e} on CG1 nodes, the slope K "
         f"was fitted against", calibration,
         np.minimum(calibration["sin_a"], cap), "nodes"),
        ("calibration half, uncapped", calibration, calibration["sin_a"],
         "nodes"),
        ("forward half, uncapped, as the forward runs today", forward,
         forward["sin_a"], "cells"),
        # A cap inside forcing.compute_sin_alpha would act on this DG0 slope,
        # so the capped forward row caps it here.
        (f"forward half, capped at {cap:.0e} on DG0 cells", forward,
         np.minimum(forward["sin_a"], cap), "cells"),
    ]

    obs_total = float(d["obs_total_gtyr"]) if "obs_total_gtyr" in d else float("nan")
    for label, g, sin_a, dofs in cases:
        floating, area = g["floating"], g["area"]
        afl = float(area[floating].sum())
        melt = np.where(floating,
                        quadratic_mixed_slope(g["tf"], g["sal"], sin_a, K=g["K"]),
                        0.0)
        over = floating & (melt > bound_m_yr)
        a_over = float(area[over].sum())
        PETSc.Sys.Print(
            f"\n  --- sin_alpha {label} ---\n"
            f"  floating                {afl / 1e6:12.1f} km^2 over "
            f"{int(floating.sum())} {dofs}\n"
            f"  melt max                {melt.max():12.1f} m/yr\n"
            f"  melt p99 (floating)     {np.quantile(melt[floating], 0.99):12.1f} m/yr\n"
            f"  melt area-mean          "
            f"{float((melt * area)[floating].sum()) / afl:12.2f} m/yr\n"
            f"  integrated              "
            f"{float((melt * area)[floating].sum()) * RHO_I / 1e12:12.0f} Gt/yr "
            f"(K was fitted to {obs_total:.0f})\n"
            f"  {dofs} past the bound    {int(over.sum()):12d}\n"
            f"  area past the bound     {a_over / 1e6:12.1f} km^2 "
            f"({100 * a_over / afl:.3f}% of floating)")

        if not over.any():
            continue
        PETSc.Sys.Print(f"  worst {dofs} (x km, y km, melt m/yr, TF K, "
                        f"draft m, sin_alpha, area km^2):")
        for i in np.argsort(-melt)[:10]:
            PETSc.Sys.Print(
                f"    {g['x'][i] / 1e3:9.1f} {g['y'][i] / 1e3:9.1f} "
                f"{melt[i]:9.1f} {g['tf'][i]:6.2f} {g['draft'][i]:8.1f} "
                f"{sin_a[i]:9.2e} {area[i] / 1e6:8.2f}")
        # A dof whose own area is a small fraction of an 8 km pixel cannot
        # fill that pixel on its own, so its value reaching the grid means the
        # pixel carried little other floating ice.
        PETSc.Sys.Print(
            f"  median area of the {dofs} past the bound: "
            f"{np.median(area[over]) / 1e6:.2f} km^2, against "
            f"{ISMIP7_DX ** 2 / 1e6:.0f} km^2 for an 8 km pixel")


if __name__ == "__main__":
    main()
