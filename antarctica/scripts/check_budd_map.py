#!/usr/bin/env python
r"""Check a Budd MAP: does any floating cell carry friction, and does the
forward reproduce the MAP's velocity?

    python antarctica/scripts/check_budd_map.py MAP.h5 [--forward]

Census (always): on the MAP's reference geometry, evaluate the production
``N_hat`` (``dual_friction.budd_nhat``, ``N_ref = N``) under the two shelf
gates and count floating cells with ``N_hat > 0``:

  old         ``conditional(N > 0, N_hat, 0)`` -- a sign test on the roundoff
              residue ``N = max(p_I - p_W, 0)``; the delta floor then lifts a
              roundoff-positive shelf cell to the cap (Sep 13 2026 finding)
  He x sign   ``He * conditional(N > 0, N_hat, 0)`` -- the first fix; leaves
              cells floating by a few metres (inside the He band) with friction
  production  ``dual_friction.budd_nhat``: ``conditional(HAF > 0, N_hat, 0)``,
              an exact gate that leaves grounded friction unscaled

The old count is the defect's footprint on that MAP; the production count
must be 0 (floating is HAF < 0, the same test the gate uses).
Also reports the grounded cells above 1 (the delta floor's footprint) and the
He band (0 < He < 1), which the He x sign row acts on.

``--forward``: run ``simulation.setup_model()`` on the MAP (the cold-start
diagnostic re-solve under the current law) and print the relative L2
distance between the solved velocity and the MAP's saved one. A MAP inverted
under the current law reproduces itself to solver tolerance (~1e-9); a MAP
inverted under the old gate does not.
"""
import argparse
import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
from icepack2_tools.dual_friction import (budd_nhat, budd_nhat_ungated,     # noqa: E402
                                          effective_pressure, grounded_mask)
from icepack2_tools.grounding import height_above_flotation               # noqa: E402
from firedrake import (CheckpointFile, Constant, Function, FunctionSpace,   # noqa: E402
                       assemble, conditional, dx, gt, inner, sqrt)
from mpi4py import MPI                                                     # noqa: E402


def gsum(comm, x):
    return comm.allreduce(float(x), op=MPI.SUM)


def census(map_path, nhat_floor, nhat_cap, gl_width):
    with CheckpointFile(map_path, "r") as chk:
        mesh = chk.load_mesh()
        H = chk.load_function(mesh, name="thickness")
        s = chk.load_function(mesh, name="surface")
        b = chk.load_function(mesh, name="bed")
        attrs = {k: chk.get_attr("/", k) for k in ("lc", "lc_coarse", "buffer_m", "mesh_basename",
                                                    "raster_sample", "misfit_norm")
                 if chk.has_attr("/", k)}
    comm = mesh.comm
    Q = FunctionSpace(mesh, "DG", 0)
    N = effective_pressure(H, s)
    He = grounded_mask(H, b, gl_width=gl_width)
    nh = budd_nhat_ungated(N, None, H, nhat_floor=nhat_floor, nhat_cap=nhat_cap)
    old = Function(Q).interpolate(conditional(gt(N, Constant(0.0)), nh, Constant(0.0)))       # sign test
    heo = Function(Q).interpolate(He * conditional(gt(N, Constant(0.0)), nh, Constant(0.0)))  # He x sign test
    new = Function(Q).interpolate(budd_nhat(N, None, H, b, nhat_floor=nhat_floor, nhat_cap=nhat_cap))
    haf = Function(Q).interpolate(height_above_flotation(H, b))
    he = Function(Q).interpolate(He)
    hq = Function(Q).interpolate(H)
    o, n, f, e, h, ho = (x.dat.data_ro for x in (old, new, haf, he, hq, heo))
    ice = h > 1.0
    floating = ice & (f < 0.0)
    grounded = ice & (f >= 0.0)
    band = ice & (e > 0.0) & (e < 1.0)
    out = {
        "cells": gsum(comm, ice.sum()),
        "floating": gsum(comm, floating.sum()),
        "floating N_hat>0 (old gate)": gsum(comm, (floating & (o > 0.0)).sum()),
        "floating N_hat at cap (old gate)": gsum(comm, (floating & (o >= nhat_cap - 1e-9)).sum()),
        "floating N_hat>0 (He x sign test)": gsum(comm, (floating & (ho > 0.0)).sum()),
        "floating N_hat>0 (production: HAF>0)": gsum(comm, (floating & (n > 0.0)).sum()),
        "grounded N_hat>1 (delta floor)": gsum(comm, (grounded & (n > 1.0 + 1e-9)).sum()),
        "grounded N_hat<1": gsum(comm, (grounded & (n < 1.0 - 1e-9)).sum()),
        "He band 0<He<1": gsum(comm, band.sum()),
    }
    return attrs, out


def node_permutation(xa, xb, decimals=6):
    r"""``(ia, ib)`` such that ``xa[ia] == xb[ib]`` row by row: the index maps
    that align two numberings of one vertex set. Raises ``SystemExit`` when
    the two are not the same set of points."""
    import numpy as np
    xa = np.asarray(xa, float)
    xb = np.asarray(xb, float)
    ia = np.lexsort(np.round(xa, decimals).T[::-1])
    ib = np.lexsort(np.round(xb, decimals).T[::-1])
    if xa.shape != xb.shape or not np.allclose(xa[ia], xb[ib], atol=10.0 ** -decimals):
        raise SystemExit("forward mesh and checkpoint mesh are not the same vertex set")
    return ia, ib


def forward_check(map_path):
    os.environ["ISMIP7_INVERSION"] = map_path
    # Force the law, do not defer to the environment: this check exists to
    # measure the Budd gate, and site_env.sh exports regularized_coulomb by
    # default, so an inherited value would re-solve the diagnostic under the
    # wrong law and return a large rel L2 for an unrelated reason.
    os.environ["ISMIP7_FRICTION"] = "budd"
    sys.path.insert(0, _ROOT)
    import simulation                                                     # noqa: E402
    ctx = simulation.setup_model()
    u = ctx["z"].subfunctions[0]
    mesh = ctx["mesh"]
    with CheckpointFile(map_path, "r") as chk:
        try:
            u_map = chk.load_function(mesh, name="velocity")
        except RuntimeError as e:
            if "velocity" in str(e):
                # the periodic (every-20-iterate) checkpoint carries the
                # controls but not the velocity; only the final save does
                raise SystemExit("no saved velocity in this MAP (a periodic checkpoint?); "
                                 "the re-solve check needs the final MAP")
            raise
        except Exception:                                                # forward mesh not from the checkpoint
            # The forward built its mesh from the .msh, the checkpoint carries
            # its own copy of the same mesh, and Firedrake numbers the two
            # differently. A raw dat copy therefore compared PERMUTED fields:
            # on a 32 km Budd MAP it reported rel L2 0.57 for a forward whose
            # nodal speeds matched the MAP's to 1e-2 m/yr once sorted (Sep 20
            # 2026). Match nodes by coordinate instead, and refuse anything
            # that is not the same vertex set.
            m2 = chk.load_mesh()
            u2 = chk.load_function(m2, name="velocity")
            if mesh.comm.size != 1:
                raise SystemExit("forward mesh is not the checkpoint mesh; run serially")
            from firedrake import SpatialCoordinate, VectorFunctionSpace
            def _coords(msh, V):
                return Function(VectorFunctionSpace(msh, V.ufl_element().family(),
                                                    V.ufl_element().degree())
                                ).interpolate(SpatialCoordinate(msh)).dat.data_ro
            ia, ib = node_permutation(_coords(mesh, u.function_space()),
                                      _coords(m2, u2.function_space()))
            u_map = Function(u.function_space())
            u_map.dat.data[ia] = u2.dat.data_ro[ib]
    num = sqrt(assemble(inner(u - u_map, u - u_map) * dx))
    den = sqrt(assemble(inner(u_map, u_map) * dx))
    return float(num / den), float(assemble(sqrt(inner(u, u)) * dx) / assemble(Constant(1.0) * dx(mesh)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("map")
    ap.add_argument("--forward", action="store_true", help="re-solve the diagnostic at the MAP geometry")
    ap.add_argument("--nhat-floor", type=float, default=float(os.environ.get("ISMIP7_BUDD_DELTA", "0.02")))
    ap.add_argument("--nhat-cap", type=float, default=float(os.environ.get("ISMIP7_BUDD_NHAT_CAP", "3.0")))
    ap.add_argument("--gl-width", type=float, default=10.0)
    a = ap.parse_args()
    rank = MPI.COMM_WORLD.rank
    attrs, out = census(a.map, a.nhat_floor, a.nhat_cap, a.gl_width)
    if rank == 0:
        print(f"MAP {os.path.basename(a.map)}  {attrs}")
        print(f"  N_hat: delta={a.nhat_floor} cap={a.nhat_cap} gl_width={a.gl_width} m (N_ref = N, the reference geometry)")
        w = max(len(k) for k in out)
        for k, v in out.items():
            print(f"  {k:<{w}}  {int(v):8d}")
    if a.forward:
        rel, umean = forward_check(a.map)
        if rank == 0:
            print(f"  forward re-solve (friction={os.environ['ISMIP7_FRICTION']}) vs MAP "
                  f"velocity: rel L2 = {rel:.3e}  (mean |u| {umean:.1f} m/yr)")


if __name__ == "__main__":
    main()
