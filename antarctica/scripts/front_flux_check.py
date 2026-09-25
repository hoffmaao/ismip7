#!/usr/bin/env python3
r"""Flux and speed at the ice front of a buffered-mesh state (issue #115).

Read-only: loads each checkpoint (a timing cache, or any state written by
``simulation.save_model_state``) and assembles facet and cell integrals, with
no solve. For each state, under its own velocity and under the
``velocity_obs`` it carries:

- the t=0 front: the facets between t=0 ice (``H_init >= 1 m``) and the t=0
  ice-free buffer, all of them and split into floating and grounded ice by
  the current height above flotation;
- the current front: the facets between cells holding more than 1 m and the
  rest, which moves if the front band empties;
- the grounding line, grounded to floating t=0 ice;
- the front band, the t=0 ice cells sharing a facet with the buffer: mass,
  mean thickness, cells emptied, and the apparent-MB reference it carries;
- the speed of floating and grounded ice binned by distance from the t=0
  front, and of thin (1 to 10 m) ice, where the floor-cell ocean drag ramps;
- the band's fluidity and friction controls against the rest of the ice.

Fluxes are upwind: the ice side's thickness times the outward normal speed,
the transport's own rule. Every number is an assembled integral or a
cross-rank reduction (AGENTS.md section 3), and no field is copied between
files, so any rank count reads any file.

    mpiexec -n 32 python antarctica/scripts/front_flux_check.py \
        --out front.json t0=CACHE.h5 step1=RUN_2500_t2015.1.h5 ...
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, os.fspath(_SCRIPTS.parents[1]))

import firedrake as fd  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402

from icepack2_tools.levelset import initial_distance  # noqa: E402
from icepack2_tools.mpi_stats import (  # noqa: E402
    global_count, global_extreme_location, global_max,
)

RHO_I = 917.0
RHO_W = 1024.0
RHO_GT = RHO_I / 1.0e12          # m^3 of ice -> Gt, as in audit_timing_cache.py
FRONT_HMIN = 1.0                  # ISMIP7_FRONT_HMIN's default
THIN_MAX = 10.0                   # ISMIP7_H_OCEAN's default: the drag ramp
DISTANCE_BINS_KM = (0.0, 5.0, 10.0, 25.0, 50.0, 100.0, 250.0, 1.0e9)
OPTIONAL = ("velocity_obs", "a_ref_mb", "log_fluidity", "log_friction")


def load(path):
    """Mesh, fields and t_yr of one checkpoint; optional fields may be None."""
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        fields = {name: chk.load_function(mesh, name=name)
                  for name in ("thickness", "H_init", "bed", "velocity")}
        for name in OPTIONAL:
            try:
                fields[name] = chk.load_function(mesh, name=name)
            except Exception:
                fields[name] = None
        t_yr = float(chk.get_attr("/", "t_yr")) if chk.has_attr("/", "t_yr") else None
    return mesh, fields, t_yr


def check(label, path):
    mesh, f, t_yr = load(path)
    comm = mesh.comm
    # DG0 whatever the stored geometry space is: the facet gates below read
    # a('+') * b('-'), which a continuous indicator makes identically zero.
    Q = fd.FunctionSpace(mesh, "DG", 0)
    h, H0, bed = (fd.Function(Q).project(f[name]) for name in ("thickness", "H_init", "bed"))
    u, uo = f["velocity"], f["velocity_obs"]

    def indicator(condition):
        return fd.Function(Q).interpolate(fd.conditional(condition, 1.0, 0.0))

    haf = h + fd.min_value(bed, 0.0) * (RHO_W / RHO_I)
    ice0 = indicator(fd.ge(H0, FRONT_HMIN))
    beyond0 = fd.Function(Q).interpolate(1.0 - ice0)
    gr0 = fd.Function(Q).interpolate(ice0 * fd.conditional(fd.gt(haf, 0.0), 1.0, 0.0))
    fl0 = fd.Function(Q).interpolate(ice0 - gr0)
    ice_now = indicator(fd.gt(h, FRONT_HMIN))
    open_now = fd.Function(Q).interpolate(1.0 - ice_now)
    fl_now = fd.Function(Q).interpolate(ice_now * fd.conditional(fd.gt(haf, 0.0), 0.0, 1.0))

    n = fd.FacetNormal(mesh)

    def over(a, b, integrand):
        """Sum over facets with an ``a`` cell on one side and a ``b`` cell on
        the other; ``integrand(side)`` is evaluated on the ``a`` side."""
        form = (a("+") * b("-") * integrand("+") + a("-") * b("+") * integrand("-")) * fd.dS
        return float(fd.assemble(form))

    def facets(a, b):
        length = over(a, b, lambda s: 1.0)
        row = {"length_km": length / 1e3}
        row["thickness_mean_m"] = over(a, b, lambda s: h(s)) / length if length else None
        for tag, v in (("model", u), ("obs", uo)):
            if v is None:
                continue
            speed = fd.sqrt(fd.dot(v, v))
            un = over(a, b, lambda s: fd.dot(v(s), n(s)))
            row[f"{tag}_mean_un_m_per_yr"] = un / length if length else None
            row[f"{tag}_mean_speed_m_per_yr"] = (
                over(a, b, lambda s: speed(s)) / length if length else None)
            row[f"{tag}_flux_gt_per_yr"] = over(
                a, b, lambda s: h(s) * fd.max_value(fd.dot(v(s), n(s)), 0.0)) * RHO_GT
        return row

    out = {"label": label, "path": os.path.abspath(path), "t_yr": t_yr,
           "mass_gt": float(fd.assemble(h * fd.dx)) * RHO_GT}
    for name, a, b in (("front_t0_all", ice0, beyond0),
                       ("front_t0_floating", fl0, beyond0),
                       ("front_t0_grounded", gr0, beyond0),
                       ("front_now_all", ice_now, open_now),
                       ("front_now_floating", fl_now, open_now),
                       ("grounding_line", gr0, fl0)):
        out[name] = facets(a, b)

    # The band: t=0 ice cells with at least one facet on the t=0 buffer.
    w = fd.TestFunction(Q)
    touch = fd.assemble((w("+") * beyond0("-") + w("-") * beyond0("+")) * fd.dS)
    band = fd.Function(Q)
    band.dat.data[:] = np.where((ice0.dat.data_ro > 0.5) & (touch.dat.data_ro > 0.0),
                                1.0, 0.0)
    band_area = float(fd.assemble(band * fd.dx))
    in_band = band.dat.data_ro > 0.5
    thin_left = h.dat.data_ro <= FRONT_HMIN
    out["band"] = {
        "cells": global_count(in_band, comm),
        "area_km2": band_area / 1e6,
        "mass_gt": float(fd.assemble(band * h * fd.dx)) * RHO_GT,
        "thickness_mean_m": float(fd.assemble(band * h * fd.dx)) / band_area,
        "cells_at_or_below_1m": global_count(in_band & thin_left, comm),
    }
    out["t0_ice_cells_at_or_below_1m"] = global_count(
        (ice0.dat.data_ro > 0.5) & thin_left, comm)
    if f["a_ref_mb"] is not None:
        a_ref = f["a_ref_mb"]
        out["band"]["a_ref_net_gt_per_yr"] = float(fd.assemble(band * a_ref * fd.dx)) * RHO_GT
        out["a_ref_net_gt_per_yr"] = float(fd.assemble(a_ref * fd.dx)) * RHO_GT

    # Speed by distance from the t=0 front, which the pinned front keeps.
    dist = initial_distance(mesh, H0, h_min=FRONT_HMIN)       # negative in ice
    depth_km = fd.Function(Q).interpolate(-dist / 1e3)
    bins = []
    for lo, hi in zip(DISTANCE_BINS_KM[:-1], DISTANCE_BINS_KM[1:]):
        within = fd.conditional(fd.And(fd.ge(depth_km, lo), fd.lt(depth_km, hi)), 1.0, 0.0)
        row = {"from_km": lo, "to_km": hi if hi < 1e8 else None}
        for region, ind in (("floating", fl0), ("grounded", gr0)):
            area = float(fd.assemble(ind * within * fd.dx))
            row[f"{region}_area_km2"] = area / 1e6
            for tag, v in (("model", u), ("obs", uo)):
                if v is None:
                    continue
                speed = fd.sqrt(fd.dot(v, v))
                row[f"{region}_{tag}_speed_m_per_yr"] = (
                    float(fd.assemble(ind * within * speed * fd.dx)) / area if area else None)
        bins.append(row)
    out["speed_by_distance"] = bins

    thin = fd.Function(Q).interpolate(
        ice0 * fd.conditional(fd.And(fd.gt(h, FRONT_HMIN), fd.lt(h, THIN_MAX)), 1.0, 0.0))
    thin_area = float(fd.assemble(thin * fd.dx))
    out["thin_ice"] = {"area_km2": thin_area / 1e6}
    for tag, v in (("model", u), ("obs", uo)):
        if v is not None and thin_area:
            speed = fd.sqrt(fd.dot(v, v))
            out["thin_ice"][f"{tag}_speed_m_per_yr"] = (
                float(fd.assemble(thin * speed * fd.dx)) / thin_area)

    # The controls the inversion set in the band, against the rest of the ice.
    rest = fd.Function(Q).interpolate(ice0 - band)
    rest_area = float(fd.assemble(rest * fd.dx))
    for name in ("log_fluidity", "log_friction"):
        if f[name] is not None:
            out[f"band_{name}_mean"] = float(fd.assemble(band * f[name] * fd.dx)) / band_area
            out[f"rest_{name}_mean"] = float(fd.assemble(rest * f[name] * fd.dx)) / rest_area

    # The fastest node on t=0 ice (a vertex of any t=0 ice cell) and the
    # fastest node of the ice-free buffer alone. Around the first, the
    # area-weighted means over the cells sharing that vertex: an ice fraction
    # of 1 is an interior node, anything less a node on the t=0 front.
    S = fd.FunctionSpace(mesh, "CG", 1)
    w_cg = fd.TestFunction(S)
    patch_area = fd.assemble(w_cg * fd.dx).dat.data_ro

    def patch_mean(field):
        g = fd.Function(S)
        g.dat.data[:] = fd.assemble(w_cg * field * fd.dx).dat.data_ro / patch_area
        return g

    speed_cg = fd.Function(S).interpolate(fd.sqrt(fd.dot(u, u)))
    out["speed_max_m_per_yr"] = global_max(speed_cg)
    ice_fraction = patch_mean(ice0)
    on_ice = ice_fraction.dat.data_ro > 0.0
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        fd.SpatialCoordinate(mesh))
    speeds = speed_cg.dat.data_ro
    for name, mask in (("ice", on_ice), ("ice_free", ~on_ice)):
        value, where = global_extreme_location(
            np.where(mask, speeds, -1.0), xy.dat.data_ro, comm=comm)
        out[f"speed_max_{name}_nodes"] = {"m_per_yr": value, "xy_m": list(where)}
    where = out["speed_max_ice_nodes"]["xy_m"]
    if where:
        for name, field in (("ice_fraction", ice0), ("thickness_m", h), ("H_init_m", H0),
                            ("bed_m", bed), ("depth_from_t0_front_km", depth_km)):
            value = patch_mean(field).at(where, dont_raise=True)
            out["speed_max_ice_nodes"][f"patch_{name}"] = (
                float(value) if value is not None else None)

    fl = out["front_t0_floating"]
    PETSc.Sys.Print(
        f"[{label}] t={t_yr} mass {out['mass_gt']:,.1f} Gt; t0 floating front "
        f"{fl['length_km']:,.0f} km, model u.n {fl.get('model_mean_un_m_per_yr', 0):+.1f} m/yr, "
        f"flux {fl.get('model_flux_gt_per_yr', 0):.1f} Gt/yr; obs u.n "
        f"{fl.get('obs_mean_un_m_per_yr') or float('nan'):+.1f}, flux "
        f"{fl.get('obs_flux_gt_per_yr') or float('nan'):.1f} Gt/yr; band "
        f"{out['band']['mass_gt']:,.1f} Gt, {out['band']['cells_at_or_below_1m']} of "
        f"{out['band']['cells']} cells at or below 1 m; GL model flux "
        f"{out['grounding_line'].get('model_flux_gt_per_yr', 0):,.1f} Gt/yr")
    return out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("states", nargs="+", help="LABEL=PATH, one per checkpoint")
    parser.add_argument("--out", required=True, help="JSON written on rank 0")
    args = parser.parse_args()
    results = []
    for spec in args.states:
        label, sep, path = spec.partition("=")
        if not sep:
            raise SystemExit(f"state {spec!r} is not LABEL=PATH")
        results.append(check(label, path))
    if fd.COMM_WORLD.rank == 0:
        with open(args.out, "w") as stream:
            json.dump(results, stream, indent=1)
    PETSc.Sys.Print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
