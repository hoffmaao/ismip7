#!/usr/bin/env python3
"""The t = 0 geometry of a mesh under each raster sampling, without a solve
(issue #167).

For each variant (``vertex``, and ``greene<year>`` with each front fill of
``geometry.FRONT_FILLS``) the BedMachine bed and thickness are put on the
mesh's DG0 cells as a cold start puts them (``geometry.sample_bed_thickness``,
then the lake rule and the surface from flotation), and a checkpoint is
written with ``thickness``, ``bed``, ``surface``, ``H_init`` and the
MEaSUREs ``velocity_obs`` for ``front_flux_check.py``, which then gives the
front band, its thickness and its flux under the observed velocity. This
prints, per variant, the front-cell counts and the ice it holds: floating and
grounded area and mass, and the cells the melt falls on.

    mpiexec -n 32 python antarctica/scripts/probe_front_cells.py \\
        --mesh antarctica/mesh/antarctica_5000_2000_buffered20000_front2015.msh \\
        --out-dir results/i167_probe [--variants vertex greene2015:empty ...]
"""
import argparse
import glob
import json
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _ROOT)

import icepack2_tools.dual_friction  # noqa: F401,E402  (icepack2 before UFL forms)

import numpy as np  # noqa: E402,F401
import firedrake as fd  # noqa: E402
from mpi4py import MPI  # noqa: E402
from firedrake import (  # noqa: E402
    Constant, Function, FunctionSpace, Mesh, TestFunction, VectorFunctionSpace,
    assemble, dx, max_value,
)
from firedrake.petsc import PETSc  # noqa: E402

from icepack2_tools.forcing import melt_receiving  # noqa: E402
from icepack2_tools.geometry import (  # noqa: E402
    FRONT_FILLS, raise_bed_to_lake_ice_base, sample_bed_thickness,
)
from icepack2_tools.runconfig import (  # noqa: E402
    lake_ice_base, obs_data_root, raster_front_year,
)

RHO_I, RHO_W = 917.0, 1024.0


def _variants(names):
    out = []
    for v in names:
        method, _, fill = v.partition(":")
        if raster_front_year(method) is None:
            out.append((method, None))
        else:
            for f in ([fill] if fill else list(FRONT_FILLS)):
                out.append((method, f))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--mesh", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--variants", nargs="+", default=["vertex", "greene2015"])
    args = p.parse_args()

    import icepack
    import rasterio

    mesh = Mesh(args.mesh, name="firedrake_default")
    comm = mesh.comm
    stem = os.path.splitext(os.path.basename(args.mesh))[0]
    if comm.rank == 0:
        os.makedirs(args.out_dir, exist_ok=True)
    Q = FunctionSpace(mesh, "CG", 1)
    Q0 = FunctionSpace(mesh, "DG", 0)
    V = VectorFunctionSpace(mesh, "CG", 1)
    root = obs_data_root()
    bm_fn = sorted(glob.glob(os.path.join(root, "bedmachine", "*.nc")))[0]
    vel_fn = sorted(glob.glob(os.path.join(root, "velocity", "*.nc")))[0]
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:VX"), rasterio.open(f"netcdf:{vel_fn}:VY")),
        V, fillvalue=0.0)
    u_obs = Function(V, name="velocity_obs").interpolate(u_obs)
    area = assemble(TestFunction(Q0) * dx).dat.data_ro
    n_cells = int(comm.allreduce(int(Q0.dof_dset.size)))
    PETSc.Sys.Print(f"{stem}: {n_cells} cells; BedMachine {os.path.basename(bm_fn)}")

    rows = []
    for method, fill in _variants(args.variants):
        label = method if fill is None else f"{method}_{fill}"
        b, H, counts = sample_bed_thickness(bm_fn, Q0, Q, floor=0.0, method=method, fill=fill)
        if lake_ice_base():
            raise_bed_to_lake_ice_base(b, H, bm_fn, Q0, Q, method=method)
        s = Function(Q0, name="surface").interpolate(
            max_value(b + H, (Constant(1.0) - Constant(RHO_I / RHO_W)) * H))
        h, bb = H.dat.data_ro, b.dat.data_ro
        ice = h >= 1.0
        afloat = (RHO_I / RHO_W) * h <= -bb
        floating = ice & afloat
        grounded = ice & ~afloat
        melt = melt_receiving(s.dat.data_ro, bb, h)

        def gsum(x):
            return comm.allreduce(float(x), op=MPI.SUM)
        row = {
            "mesh": stem, "variant": label, "cells": n_cells, "counts": counts,
            "floating_km2": round(gsum(area[floating].sum()) / 1e6, 1),
            "floating_gt": round(gsum((area * h)[floating].sum()) * RHO_I / 1e12, 1),
            "grounded_km2": round(gsum(area[grounded].sum()) / 1e6, 1),
            "grounded_gt": round(gsum((area * h)[grounded].sum()) * RHO_I / 1e12, 1),
            "melt_cells": int(gsum(melt.sum())),
            "melt_km2": round(gsum(area[melt].sum()) / 1e6, 1),
        }
        rows.append(row)
        PETSc.Sys.Print(json.dumps(row))
        H.rename("thickness")
        b.rename("bed")
        H_init = Function(Q0, name="H_init").assign(H)
        out = os.path.join(args.out_dir, f"{stem}_{label}.h5")
        with fd.CheckpointFile(out, "w") as chk:
            chk.save_mesh(mesh)
            for f in (H, b, s, H_init, u_obs):
                chk.save_function(f)
            chk.set_attr("/", "raster_sample", method)
            chk.set_attr("/", "front_fill", fill or "none")
            chk.set_attr("/", "mesh_basename", stem)
        PETSc.Sys.Print(f"  wrote {out}")
    if comm.rank == 0:
        with open(os.path.join(args.out_dir, f"{stem}_probe.json"), "w") as f:
            json.dump(rows, f, indent=1)


if __name__ == "__main__":
    main()
