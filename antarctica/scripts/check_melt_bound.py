#!/usr/bin/env python
r"""How far does the draft-slope cap move melt past the variable request's bound?

    python antarctica/scripts/check_melt_bound.py [--npz antarctica/results/calibrated_K_per_basin_2000.npz]
        [--ocx [main|cold|warm|vary]] [--ocx-years 2000,2015,2025] [--ocx-tol 0.25]

The ISMIP7 variable request gives ``libmassbffl`` an AIS minimum of
-0.008 kg m-2 s-1 with severity ``error``. In ice-equivalent thickness that is
275.3 m/yr, and a 10-year adaptive-mesh ssp585 reached -0.0117 (402.6 m/yr) on
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

``--ocx`` adds the check core 11 needs before it runs on the ISMIP7 OCX
product. Every K here is fitted to the OI climatology, and the OCX ocean is a
different field: discussion #48 (17 September 2026) reports the OCX ``main``
thermal forcing so far from the Zhou climatology around Mertz that the
region's melt halves against a calibration made on the climatology, possibly
because OCX was built from an older extrapolated climatology. So the forward
half is melted a second time, uncapped as the forward runs, with thermal
forcing and salinity from the OCX ocean through the forward's own reader, for
each of ``--ocx-years``, and the two melts are set side by side per IMBIE2
basin and per 256 km block, the blocks because a basin total dilutes one
shelf. A basin or block whose OCX melt is off the climatology's by more than
``--ocx-tol`` is flagged and the exit status is 1, so a chain can stop on it.

Serial. Reuses calibrate_melt's loaders, so it needs the same inputs: a MAP for
the mesh, the OI climatology, the IMBIE2 basins and BedMachine.

Measured on the adaptive 2 km mesh, September 2026, with
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

With calibrated_K_per_basin_2500.npz, the coefficient file the 10-year ssp585
run read, everything else unchanged:

* calibration half, capped: maximum 54.0 m/yr, area mean 0.61, 841 Gt/yr, and
  nothing past the bound;
* calibration half, uncapped: maximum 1522.7 m/yr, area mean 3.34, 4634 Gt/yr,
  and 320 nodes past the bound over 2193.7 km2;
* forward half, uncapped: maximum 869.8 m/yr, area mean 0.92, 1380 Gt/yr, and
  63 cells past the bound over 248.6 km2 with a median cell area of 3.54 km2;
* forward half, capped: maximum 43.9 m/yr, area mean 0.34, 510 Gt/yr, and
  nothing past the bound.

The calibration half reproduces the target its own K was fitted to, 841
against 865 Gt/yr for the 2500 file, which is the internal consistency check.
The ratio of forward to calibration is the durable result, stable across both
coefficient files: 1.62 with the 2000 file and 1.64 with the 2500 file. The
10-year run booked 1860 Gt/yr with the 2500 file, above the 1380 Gt/yr its
forward half applies at the reference state, since that run carries warmer
ssp585 thermal forcing over evolving geometry; the gap is a consistent
residual. An earlier claim that the forward half lands within 7% of that run
compared two different coefficient files and does not hold.

Capping the forward's slope undershoots the target by 39% with the 2000 file,
so neither convention on its own reconciles the halves, which also differ in
floating area, mask and quadrature. Recalibrating K through the forward's own
cell by cell melt path, under whichever slope convention is chosen, reconciles
them by construction.

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
                                    compute_sin_alpha, ISMIP7Ocean,
                                    OCX, OCX_OCEAN_VARIANTS,
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


# Side of the blocks the OCX comparison also sums over: 32 pixels of the 8 km
# forcing grid, a few shelves wide, so one shelf is not lost in a basin total.
BLOCK_M = 256.0e3


def block_ids(x, y, size=BLOCK_M):
    r"""An integer id per point naming the ``size`` square it falls in, and
    the ``{id: (x centre, y centre)}`` of the squares that occur."""
    ix, iy = np.floor(x / size).astype(int), np.floor(y / size).astype(int)
    ids = ix * 100000 + iy
    centres = {int(i): ((a + 0.5) * size, (b + 0.5) * size)
               for i, a, b in zip(ids, ix, iy)}
    return ids, centres


def melt_by_group(groups, melt, area, floating):
    r"""``{group: Gt/yr}`` of ``melt`` (m/yr of ice) over the floating dofs."""
    gt = np.where(floating, melt * area, 0.0) * RHO_I / 1e12
    return {int(g): float(gt[groups == g].sum()) for g in np.unique(groups[floating])}


def off_by_more_than(reference, other, tol, floor_gt):
    r"""The groups whose ``other`` melt is off ``reference`` by more than
    ``tol`` (a fraction). Groups melting less than ``floor_gt`` either way are
    left out: a ratio of two near-zero totals flags nothing worth reading."""
    flagged = []
    for g, ref in reference.items():
        new = other.get(g, 0.0)
        if max(ref, new) < floor_gt:
            continue
        if ref <= 0.0 or abs(new / ref - 1.0) > tol:
            flagged.append(g)
    return flagged


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", default=None,
                    help="calibrated_K_per_basin_<lc>.npz (default: the one "
                         "beside this run's lc)")
    ap.add_argument("--ocx", nargs="?", const="main", default=None,
                    choices=OCX_OCEAN_VARIANTS,
                    help="also melt the forward half with this OCX ocean and "
                         "compare it with the climatology K was fitted to "
                         "(discussion #48); exit 1 if any region is off")
    ap.add_argument("--ocx-years", default="2000,2015,2025",
                    help="comma list of OCX years to compare (default %(default)s)")
    ap.add_argument("--ocx-tol", type=float, default=0.25,
                    help="flag a region whose OCX melt is off by more than "
                         "this fraction (default %(default)s)")
    ap.add_argument("--ocx-floor", type=float, default=2.0,
                    help="ignore regions melting less than this many Gt/yr "
                         "either way (default %(default)s)")
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
    if "obs_total_gtyr" in d:
        target = f" (K was fitted to {float(d['obs_total_gtyr']):.0f})"
    else:
        target = ""
        PETSc.Sys.Print(f"  {os.path.basename(npz_path)} records no "
                        f"obs_total_gtyr, so the rows print no target")

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
            f"{float((melt * area)[floating].sum()) * RHO_I / 1e12:12.0f} Gt/yr"
            f"{target}\n"
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

    if a.ocx is None:
        return 0
    return compare_with_ocx(a, forward)


def compare_with_ocx(a, g):
    r"""Melt the forward half with the OCX ocean and set it beside the melt
    from the climatology K was fitted to. Returns the exit status."""
    ocean = ISMIP7Ocean(scenario=OCX, variant=a.ocx)
    years = [int(v) for v in a.ocx_years.split(",") if v.strip()]
    floating, area, sin_a = g["floating"], g["area"], g["sin_a"]
    basin = np.round(cm._grid_interp(cm.IMBIE2_NC, "basinNumber", g["x"], g["y"])).astype(int)
    block, centres = block_ids(g["x"], g["y"])

    def melt_with(tf, sal):
        return np.where(floating, quadratic_mixed_slope(tf, sal, sin_a, K=g["K"]), 0.0)

    reference = melt_with(g["tf"], g["sal"])
    PETSc.Sys.Print(
        f"\n=== OCX ocean '{a.ocx}' against the climatology K was fitted to "
        f"(discussion #48) ===\n"
        f"  forward half, uncapped; a region is flagged past "
        f"{100 * a.ocx_tol:.0f}%, regions under {a.ocx_floor:g} Gt/yr ignored")
    status = 0
    for year in years:
        other = melt_with(ocean.get_thermal_forcing(year, g["x"], g["y"], draft=g["draft"]),
                          ocean.get_salinity(year, g["x"], g["y"], draft=g["draft"]))
        for name, groups, where in (("IMBIE2 basin", basin, None),
                                    (f"{BLOCK_M / 1e3:.0f} km block", block, centres)):
            ref = melt_by_group(groups, reference, area, floating)
            new = melt_by_group(groups, other, area, floating)
            bad = off_by_more_than(ref, new, a.ocx_tol, a.ocx_floor)
            PETSc.Sys.Print(
                f"\n  --- {year}, by {name}: climatology {sum(ref.values()):.0f} Gt/yr, "
                f"OCX {sum(new.values()):.0f} Gt/yr, {len(bad)} of {len(ref)} flagged ---")
            # every basin, but only the blocks that are off: there are hundreds
            for grp in (sorted(ref) if where is None else
                        sorted(bad, key=lambda b: -abs(new.get(b, 0.0) - ref[b]))):
                at = (f"{grp:6d}" if where is None else
                      f"x {where[grp][0] / 1e3:7.0f} km  y {where[grp][1] / 1e3:7.0f} km")
                ratio = new.get(grp, 0.0) / ref[grp] if ref[grp] > 0 else float("inf")
                PETSc.Sys.Print(
                    f"    {at}  climatology {ref[grp]:8.1f}  OCX {new.get(grp, 0.0):8.1f} Gt/yr"
                    f"  ratio {ratio:5.2f}{'   <-- OFF' if grp in bad else ''}")
            if bad:
                status = 1
    if status:
        PETSc.Sys.Print(
            "\n  The OCX ocean and the climatology disagree by more than the "
            "tolerance somewhere.\n  K is fitted to the climatology, so a "
            "protocol-forced core 11 melts those regions\n  differently from "
            "its own calibration. See discussion #48 before running it.")
    return status


if __name__ == "__main__":
    sys.exit(main())
