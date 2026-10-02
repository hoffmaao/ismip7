#!/usr/bin/env python3
r"""Cost of the inversion's per-evaluation prior work at production sizes.

Issue 156 measured about 27 s of each 1 km L-BFGS-B evaluation outside the
forward and the adjoint (13 s at 2 km), the same under full_mumps and
scpc_gamg and the same on 16, 32 and 64 ranks. The spans
inversion_icepack2.py now records put that time on the 32 km mesh at 0.12 s,
too small to say how it grows, so this refines the same mesh uniformly
(``MeshHierarchy``, four times the vertices a level; level 4 is 1.6 million,
the 1 km mesh has 1.87 million) and times, per level, the work an evaluation
does under ``ISMIP7_PRIOR_FORM=bilaplacian`` that is neither solve.

Variants, each run per control (theta, then phi) as an evaluation runs it:

``newton``
    the driver's code: ``clear_caches()``, then ``EquationSolver(M f - A ctrl
    == 0, f)`` with no solver parameters, so Firedrake's defaults: a Newton
    solve, MUMPS LU, and MUMPS's default sequential analysis (ICNTL(28)=1).
``newton_pt``
    the same with MUMPS parallel analysis (ICNTL(28)=2, PT-Scotch), which
    the state solves' MUMPS options already set.
``factored``
    ``M`` factored once (MUMPS Cholesky, the options of the driver's
    ``_prior_metric_solvers``); an evaluation assembles ``A ctrl`` and
    back-substitutes. The one-time factorisation is reported apart.
``cg``
    CG with Jacobi on ``M`` to rtol 1e-12, no factorisation.

Each variant also assembles the energy and the gradient ``A f``. The
replicated gathers of ``objective_and_gradient`` (``func_to_global`` four
times, ``global_to_func`` twice) are timed as well. Every number is the
slowest rank's (``icepack2_tools.profiling.Spans``).

usage:
    mpiexec -n P python antarctica/scripts/probe_eval_overhead.py MESH.msh \
        [--levels 4] [--repeats 3] [--json out.json]
"""

import argparse
import json
import os
import sys

import numpy as np

import firedrake as fd
from firedrake import (COMM_WORLD, Function, FunctionSpace, LinearSolver,
                       MeshHierarchy, SpatialCoordinate, TestFunction,
                       TrialFunction, assemble, cos, dx, inner, sin)
from firedrake.petsc import PETSc
from tlm_adjoint.firedrake import EquationSolver, clear_caches

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from icepack2_tools.prior import (  # noqa: E402
    bilaplacian_aux_residual, bilaplacian_coeffs, bilaplacian_energy_form,
    prior_operator_form)
from icepack2_tools.profiling import Spans  # noqa: E402

# The driver's settings: fc_params, and the bi-Laplacian prior of the issue
# 156 runs (sigma 0.3, rho 7500 m for both controls).
FC = {"quadrature_degree": 4}
DELTA, GAMMA = bilaplacian_coeffs(0.3, 7500.0)
CHOLESKY = {"ksp_type": "preonly", "pc_type": "cholesky",
            "pc_factor_mat_solver_type": "mumps"}
NEWTON_PT = {"snes_type": "newtonls", "snes_linesearch_type": "basic",
             "ksp_type": "preonly", "pc_type": "lu",
             "pc_factor_mat_solver_type": "mumps",
             "mat_mumps_icntl_28": 2, "mat_mumps_icntl_29": 1}
CG = {"ksp_type": "cg", "pc_type": "jacobi", "ksp_rtol": 1e-12,
      "ksp_atol": 0.0}
VARIANTS = ("newton", "newton_pt", "factored", "cg")


# Copies of the driver's MPI helpers (closures in inversion_icepack2.main).
def func_to_global(f):
    with f.dat.vec_ro as v:
        scatter, x_seq = PETSc.Scatter.toAll(v)
        scatter.scatter(v, x_seq, mode=PETSc.Scatter.Mode.FORWARD)
        result = x_seq.array.copy()
        scatter.destroy()
        x_seq.destroy()
        return result


def global_to_func(arr, f):
    with f.dat.vec_wo as v:
        x_seq = PETSc.Vec().createSeq(len(arr), comm=PETSc.COMM_SELF)
        x_seq.array[:] = arr
        scatter, _ = PETSc.Scatter.toAll(v)
        scatter.scatter(x_seq, v, mode=PETSc.Scatter.Mode.REVERSE)
        scatter.destroy()
        x_seq.destroy()


def probe_level(mesh, repeats):
    Q = FunctionSpace(mesh, "CG", 1)
    test, trial = TestFunction(Q), TrialFunction(Q)
    x, y = SpatialCoordinate(mesh)
    shape = {"theta": 0.3 * sin(x / 2.0e5) * cos(y / 3.0e5),
             "phi": 0.3 * cos(x / 2.5e5) * sin(y / 1.5e5)}
    ctrl = {k: Function(Q, name=k) for k in shape}
    aux = {v: {k: Function(Q) for k in shape} for v in VARIANTS}
    spans = Spans(COMM_WORLD)

    # The two persistent solvers, built and factored once (timed apart).
    M = assemble(inner(trial, test) * dx, form_compiler_parameters=FC)
    solvers = {"factored": LinearSolver(M, solver_parameters=CHOLESKY),
               "cg": LinearSolver(M, solver_parameters=CG)}
    setup = {}
    for name, slv in solvers.items():
        f = Function(Q)
        with spans(name):
            slv.solve(f, assemble(inner(Function(Q).assign(1.0), test) * dx))
        setup.update(spans.reduce())

    rows = []
    for r in range(repeats):
        # A new control vector each repeat, as each evaluation brings one.
        for k, expr in shape.items():
            ctrl[k].interpolate((1.0 + 0.05 * r) * expr)
        row = {"repeat": r, "energy": {}, "cg_iterations": []}
        for v in VARIANTS:
            if v.startswith("newton"):
                clear_caches()  # what forward() does at every evaluation
            energy = 0.0
            for k in ("theta", "phi"):
                f = aux[v][k]
                with spans(f"{v}.solve"):
                    if v == "newton":
                        EquationSolver(
                            bilaplacian_aux_residual(ctrl[k], f, test, DELTA, GAMMA) == 0,
                            f, form_compiler_parameters=FC).solve()
                    elif v == "newton_pt":
                        EquationSolver(
                            bilaplacian_aux_residual(ctrl[k], f, test, DELTA, GAMMA) == 0,
                            f, form_compiler_parameters=FC,
                            solver_parameters=NEWTON_PT).solve()
                    else:
                        b = assemble(prior_operator_form(ctrl[k], test, DELTA, GAMMA),
                                     form_compiler_parameters=FC)
                        solvers[v].solve(f, b)
                        if v == "cg":
                            row["cg_iterations"].append(
                                solvers[v].ksp.getIterationNumber())
                with spans(f"{v}.energy"):
                    energy += float(assemble(bilaplacian_energy_form(f),
                                             form_compiler_parameters=FC))
                with spans(f"{v}.grad"):
                    assemble(prior_operator_form(f, test, DELTA, GAMMA),
                             form_compiler_parameters=FC)
            row["energy"][v] = energy
        scratch = Function(Q)
        with spans("gather.func_to_global_x4"):
            parts = [func_to_global(ctrl["theta"]), func_to_global(ctrl["phi"]),
                     func_to_global(ctrl["theta"]), func_to_global(ctrl["phi"])]
        with spans("gather.global_to_func_x2"):
            global_to_func(parts[0], scratch)
            global_to_func(parts[1], scratch)
        row["seconds"] = spans.reduce()
        rows.append(row)
    return {"vertices": int(Q.dim()), "setup_seconds": setup, "repeats": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mesh")
    ap.add_argument("--levels", type=int, default=4)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--json", default="")
    args = ap.parse_args()

    base = fd.Mesh(args.mesh)
    hierarchy = MeshHierarchy(base, args.levels)
    out = {"mesh": os.path.basename(args.mesh), "ranks": COMM_WORLD.size,
           "levels": []}
    for level, mesh in enumerate(hierarchy):
        res = probe_level(mesh, args.repeats)
        res["level"] = level
        out["levels"].append(res)
        # the first repeat carries form compilation; report the later ones
        later = res["repeats"][1:] or res["repeats"]
        med = {k: float(np.median([r["seconds"][k] for r in later]))
               for k in later[0]["seconds"]}
        e = later[-1]["energy"]
        PETSc.Sys.Print(
            f"level {level}: {res['vertices']} vertices on {COMM_WORLD.size} ranks; "
            f"setup factored={res['setup_seconds']['factored']:.3f}s "
            f"cg={res['setup_seconds']['cg']:.3f}s; "
            f"cg iterations {later[-1]['cg_iterations']}")
        for v in VARIANTS:
            parts = {s: med[f"{v}.{s}"] for s in ("solve", "energy", "grad")}
            PETSc.Sys.Print(
                f"  {v:10s} {sum(parts.values()):8.3f}s  "
                + " ".join(f"{s}={t:.3f}" for s, t in parts.items())
                + f"  energy rel. diff from newton "
                f"{abs(e[v] - e['newton']) / abs(e['newton']):.1e}")
        PETSc.Sys.Print(
            "  gathers    "
            + " ".join(f"{k.split('.', 1)[1]}={t:.3f}s"
                       for k, t in med.items() if k.startswith("gather.")))
    if args.json and COMM_WORLD.rank == 0:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
