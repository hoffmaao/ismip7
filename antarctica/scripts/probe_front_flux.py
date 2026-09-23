#!/usr/bin/env python3
r"""Where the ice front sits, and how much ice crosses it.

    python antarctica/scripts/probe_front_flux.py STATE.h5
    python antarctica/scripts/probe_front_flux.py STATE.h5 --hmin 1,50,100,150

A calving law can only remove ice the front actually runs through, so the
threshold that defines the front, ``ISMIP7_FRONT_HMIN``, decides which ice the
law sees. That threshold is one metre by default, and BedMachine averaged onto
a cell smears the coastline, so the outermost ice cells hold a fraction of a
calving face's thickness. A law tuned for a 200 m front and applied to a
40 m fringe removes almost nothing whatever its threshold, and the run's own
calving tally and the submitted ``licalvf`` read near zero with it.

This reports, per threshold, the front's extent, its length-weighted
thickness, the outward normal speed and the flux across it, so the choice of
threshold is made against numbers rather than left at its default. It reads a
checkpoint and solves nothing.

The outward flux is ``sum_front max(u . n, 0) h L``, with ``n`` the level
set's own unit gradient and ``L`` the front length in the cell. That is the
ice the level set can act on, and it counts every facet between ice and
ice-free cells, so a thin patch inside the sheet that ice flows through
contributes its inflow side. The signed ``sum_front (u . n) h L`` cancels
that through-flux, which makes the net the column to compare with the
observed calving flux, and the gap between the two measures how much of the
front is interior.

``ISMIP7_FRONT_HMIN`` also sets the t=0 extent mask, the sliver removal
threshold and the open-water classification in ``simulation.py``, so a
threshold chosen here changes the run's t=0 mass and books every cell that
thins below it as calving. Audit a run at a raised threshold (mass budget,
``check_ismip6_track.py``, t=0 mass) before adopting it.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from firedrake import CheckpointFile, Function, FunctionSpace, dot
from firedrake.petsc import PETSc

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

from icepack2_tools.levelset import LevelSet
from icepack2_tools.mpi_stats import global_range

#: kg of ice per m^3, over 1e12, so a volume in m^3 becomes Gt.
RHO_GT = 917.0 / 1e12

#: Observed Antarctic calving flux, for the one comparison that matters
#: (Rignot et al. 2013, 1265 +/- 140 Gt/yr).
OBSERVED_CALVING_GT_YR = 1265.0

DEFAULT_THRESHOLDS = (1.0, 10.0, 25.0, 50.0, 100.0, 150.0)


def front_flux(mesh, h, u, h_min):
    r"""``(cells, length_m, h_front_m, u_normal_m_yr, flux_gt_yr,
    net_gt_yr)`` for the front at ``h_min``, all reduced over the
    communicator: ``flux`` sums the outward part of ``(u . n) h L`` and
    ``net`` the signed total."""
    Q0 = FunctionSpace(mesh, "DG", 0)
    ls = LevelSet(mesh, h, law="none", h_min=float(h_min))
    ls._update_unit_gradient()
    un = Function(Q0).interpolate(dot(u, ls.ghat)).dat.data_ro
    length = ls.front_len.dat.data_ro
    thickness = h.dat.data_ro
    front = length > 0.0

    comm = mesh.comm

    def total(values):
        return comm.allreduce(float(np.sum(values)))

    cells = int(comm.allreduce(int(front.sum())))
    total_length = total(length[front])
    if cells == 0 or total_length == 0.0:
        return cells, 0.0, 0.0, 0.0, 0.0, 0.0
    # length-weighted, so a long thin stretch counts for what it carries
    h_front = total(thickness[front] * length[front]) / total_length
    u_front = total(un[front] * length[front]) / total_length
    flux = total((np.maximum(un, 0.0) * thickness * length)[front]) * RHO_GT
    net = total((un * thickness * length)[front]) * RHO_GT
    return cells, total_length, h_front, u_front, flux, net


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("state", help="a checkpoint carrying thickness and velocity")
    p.add_argument("--hmin", default=",".join(f"{v:g}" for v in DEFAULT_THRESHOLDS),
                   help="comma list of front thresholds in metres")
    a = p.parse_args(argv)

    with CheckpointFile(a.state, "r") as chk:
        mesh = chk.load_mesh()
        try:
            h = chk.load_function(mesh, name="thickness")
            u = chk.load_function(mesh, name="velocity")
        except RuntimeError as e:
            raise SystemExit(
                f"{a.state}: {' '.join(str(e).split())}\na periodic "
                f"inversion checkpoint has no velocity, so use a forward "
                f"state or a MAP's final save") from e
        t_yr = (float(chk.get_attr("/", "t_yr"))
                if chk.has_attr("/", "t_yr") else None)
        mesh_name = (str(chk.get_attr("/", "mesh_basename"))
                     if chk.has_attr("/", "mesh_basename") else "")

    if h.function_space().ufl_element().family() != "Discontinuous Lagrange":
        raise SystemExit(
            "the front is a cell-wise object; this needs the DG0 geometry "
            "(ISMIP7_GEOMETRY_SPACE=dg0)")

    lo, hi = global_range(h)
    PETSc.Sys.Print(f"state {os.path.basename(a.state)}"
                    + (f"  t_yr={t_yr:g}" if t_yr is not None else "")
                    + (f"  mesh {mesh_name}" if mesh_name else ""))
    PETSc.Sys.Print(f"  thickness [{lo:.0f}, {hi:.0f}] m on "
                    f"{mesh.comm.allreduce(h.dat.data_ro.size)} cells")
    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"  {'h_min':>7s} {'cells':>8s} {'front km':>9s} "
                    f"{'h_front':>8s} {'u_n':>7s} {'outward':>9s} {'net':>9s} "
                    f"{'vs obs':>7s}")
    PETSc.Sys.Print(f"  {'m':>7s} {'':>8s} {'':>9s} {'m':>8s} "
                    f"{'m/yr':>7s} {'Gt/yr':>9s} {'Gt/yr':>9s} {'':>7s}")
    for value in [float(v) for v in a.hmin.split(",") if v.strip()]:
        cells, length, h_front, u_front, flux, net = front_flux(
            mesh, h, u, value)
        if cells == 0:
            PETSc.Sys.Print(f"  {value:7.0f}  no interior front at this "
                            f"threshold")
            continue
        PETSc.Sys.Print(
            f"  {value:7.0f} {cells:8d} {length / 1e3:9,.0f} {h_front:8.0f} "
            f"{u_front:7.0f} {flux:9,.0f} {net:9,.0f} "
            f"{net / OBSERVED_CALVING_GT_YR:6.2f}x")
    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"  'vs obs' is the net flux over "
                    f"{OBSERVED_CALVING_GT_YR:.0f} Gt/yr, the observed "
                    f"Antarctic calving flux (Rignot et al. 2013).")
    PETSc.Sys.Print("  'outward' well above 'net' means much of the front "
                    "encloses interior thin patches.")
    PETSc.Sys.Print("  ISMIP7_FRONT_HMIN also sets the t=0 extent and the "
                    "sliver removal; audit a run at a")
    PETSc.Sys.Print("  raised threshold (mass budget, check_ismip6_track.py, "
                    "t=0 mass) before adopting it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
