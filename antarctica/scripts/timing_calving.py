#!/usr/bin/env python3
r"""Wall time of the calving front solver on a forward state, law by law.

A forward advances the level-set front at every transport advance, so under
adaptive substeps (4 to 64 a step) the front solver runs that many times a
step. This times what one advance does, on the state of a real forward
checkpoint and the mesh it was written on (a 2 km buffered state is the
case it is built for):

    construct   building the level set and its solvers (once per run)
    eikonal     the exact distance to the front facets of the current extent
                (every front segment gathered to every rank, then a nearest
                segment search per cell)
    rate        the law's rate interpolated on the ice cells
    advance     ``LevelSet.advance``: the two above plus the normal-flow solve
    masks       ``calving_masks`` and ``beyond_front``, which the forward reads
                after every advance

each the slowest rank's time, over ``--reps`` repetitions on the same state
(the extent does not change between them, so every repetition does the same
work). The level set, the law, the front state and its rate are built
exactly as ``simulation.run_simulation`` builds them, from the same knobs:
the geometry space, ``ISMIP7_FRONT_HMIN``, the drag gate of the state.

    mpiexec -n 32 python antarctica/scripts/timing_calving.py STATE.h5 \
        --laws "hfb:sigma_max=0.15" "vonmises:sigma_max_fl=0.275,sigma_max_gr=1.0" \
               "eigencalving:K=1.1e8" nocalving fixed \
        --reps 5 --out timing_calving.json

STATE.h5 is any checkpoint a forward wrote (ISMIP7_RESTART semantics: the
run's environment, mesh and MAP knobs as the forward that wrote it), or
``cold`` for the state a forward cold-starts from: the MAP that
ISMIP7_INVERSION names on the mesh ISMIP7_MESH names, as setup_model builds
it (on another mesh, the transfer and the startup continuation).
"""
import argparse
import json
import os
import resource
import sys
from time import perf_counter

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_SCRIPTS))
sys.path.insert(0, _SCRIPTS)

import numpy as np  # noqa: E402
from mpi4py import MPI  # noqa: E402

import firedrake as fd  # noqa: E402
from firedrake import Function, dx  # noqa: E402

from simulation import setup_model, calving_front_state, PETSc  # noqa: E402
from icepack2_tools.runconfig import front_hmin as _front_hmin  # noqa: E402


def _parse_law(spec):
    name, _, params = spec.partition(":")
    out = {}
    for kv in filter(None, params.split(",")):
        k, _, v = kv.partition("=")
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = {"true": True, "false": False}.get(v.lower(), v)
    return name, out


def _slowest(comm, seconds):
    return comm.allreduce(float(seconds), op=MPI.MAX)


def _rss_gib(comm):
    kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(kib / 1024.0 ** 2, op=MPI.MAX)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("state", help="a forward checkpoint, or cold")
    ap.add_argument("--laws", nargs="+",
                    default=["fixed", "hfb:sigma_max=0.15",
                             "vonmises:sigma_max_fl=0.275,sigma_max_gr=1.0",
                             "eigencalving:K=1.1e8", "nocalving"])
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--dt", type=float, default=0.05 / 4,
                    help="the advance's step [yr]: a 0.05 yr step at 4 substeps by default")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    comm = MPI.COMM_WORLD

    from icepack_tools import calving
    from icepack2_tools.levelset import LevelSet

    t0 = perf_counter()
    cold = a.state == "cold"
    ctx = setup_model(restart_from=None if cold else a.state)
    t_setup = _slowest(comm, perf_counter() - t0)
    mesh = ctx["mesh"]
    geom_dg1 = bool(ctx.get("geom_dg1", False))
    z = ctx["z"]
    u = z.subfunctions[0]
    if geom_dg1:
        cell_view = ctx["cell_view"]
        h_front = cell_view.mean_function(ctx["h"])
        b_front = cell_view.mean_function(ctx["b"])
    else:
        h_front, b_front = ctx["h"], ctx["b"]
    hmin = _front_hmin()
    # owned entities, summed: cells from the DG0 space, vertices from the CG1 one
    n_cells = int(comm.allreduce(int(fd.FunctionSpace(mesh, "DG", 0).dof_dset.size)))
    n_vertices = int(comm.allreduce(int(fd.FunctionSpace(mesh, "CG", 1).dof_dset.size)))
    t_state = 2015.0 if cold else float(ctx.get("t_restart") or 0.0)
    PETSc.Sys.Print(
        f"calving solver timing: {'cold start from ' + os.path.basename(os.environ.get('ISMIP7_INVERSION', '?')) if cold else os.path.basename(a.state)} at t={t_state:g}, "
        f"{n_vertices} vertices, {n_cells} cells, {comm.size} ranks, "
        f"geometry {'DG1 (front on cell means)' if geom_dg1 else 'DG0'}, "
        f"h_min {hmin:g} m, dt {a.dt:g} yr, {a.reps} repetitions; setup {t_setup:.1f} s")

    rows = []
    for spec in a.laws:
        name, params = _parse_law(spec)
        law = calving.make(name, **params)
        comm.Barrier()
        t0 = perf_counter()
        ls = LevelSet(mesh, h_front, law=law.front_mode, h_min=hmin,
                      drag_mask=ctx.get("drag_mask"))
        t_construct = _slowest(comm, perf_counter() - t0)
        rate = None
        if law.front_mode == "prescribed":
            state = calving_front_state(z, h_front, b_front, ls, ctx.get("A_map"),
                                        ctx["n_flow_val"])
            rate = law.rate(state, t_state)
        t_eik, t_rate, t_adv, t_mask = [], [], [], []
        n_front, c_front, n_beyond = 0, 0.0, 0
        for _ in range(a.reps):
            comm.Barrier()
            t0 = perf_counter()
            n_front = ls.solve_eikonal_from_extent()
            t_eik.append(_slowest(comm, perf_counter() - t0))
            if rate is not None:
                comm.Barrier()
                t0 = perf_counter()
                ls.c_cell.interpolate(ls.chi * ls._rate_expr(rate))
                t_rate.append(_slowest(comm, perf_counter() - t0))
            comm.Barrier()
            t0 = perf_counter()
            c_front = ls.advance(a.dt, u, rate=rate)
            t_adv.append(_slowest(comm, perf_counter() - t0))
            comm.Barrier()
            t0 = perf_counter()
            lsb, _frac = ls.calving_masks()
            beyond = ls.beyond_front()
            t_mask.append(_slowest(comm, perf_counter() - t0))
            n_beyond = comm.allreduce(int(np.count_nonzero(lsb)))
        row = {
            "law": spec, "front_mode": law.front_mode, "front_facets": int(n_front),
            "construct_s": t_construct,
            "eikonal_s": float(np.median(t_eik)),
            "rate_s": float(np.median(t_rate)) if t_rate else 0.0,
            "advance_s": float(np.median(t_adv)), "advance_max_s": float(np.max(t_adv)),
            "masks_s": float(np.median(t_mask)),
            "per_advance_s": float(np.median(t_adv)) + float(np.median(t_mask)),
            "mean_front_rate_m_yr": float(c_front), "cells_passed": int(n_beyond),
            "rss_gib": _rss_gib(comm),
        }
        rows.append(row)
        PETSc.Sys.Print(
            f"  {spec:46s} construct {t_construct:6.2f} s | per advance "
            f"{row['per_advance_s']:6.3f} s (eikonal {row['eikonal_s']:.3f}, "
            f"rate {row['rate_s']:.3f}, advance {row['advance_s']:.3f}, "
            f"masks {row['masks_s']:.3f}) | {n_front} front facets, mean rate "
            f"{c_front:.1f} m/yr, {n_beyond} cells passed | RSS {row['rss_gib']:.2f} GiB")
    if a.out and comm.rank == 0:
        with open(a.out, "w") as fh:
            json.dump({"state": "cold" if cold else os.path.realpath(a.state), "t": t_state,
                       "map": os.environ.get("ISMIP7_INVERSION"),
                       "mesh": os.environ.get("ISMIP7_MESH"),
                       "ranks": comm.size, "vertices": n_vertices, "cells": n_cells,
                       "geometry_space": "dg1" if geom_dg1 else "dg0", "h_min": hmin,
                       "dt": a.dt, "reps": a.reps, "setup_s": t_setup, "laws": rows},
                      fh, indent=1)
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
