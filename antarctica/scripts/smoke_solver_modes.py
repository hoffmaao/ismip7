#!/usr/bin/env python3
r"""Tiny distributed initialization test for transient solver modes.

This is not a physics or convergence qualification.  It exercises the actual
CG1-vector / DG0-symmetric-tensor / DG0-vector mixed layout, exact local
condensation, GAMG/MUMPS setup, and persistent transport-solver reuse without
loading Antarctic data.
"""

import argparse
import os
import sys

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)

import firedrake as fd  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402

from icepack2_tools.preconditioners import frozen_linearization  # noqa: E402
from icepack2_tools.solverconfig import (  # noqa: E402
    diagnostic_solver_parameters,
    linearization_state,
    transport_solver_parameters,
)
from icepack2_tools.mpi_stats import global_extreme_location  # noqa: E402


def mixed_problem(mesh, cubic_drag=0.0):
    velocity = fd.VectorFunctionSpace(mesh, "CG", 1)
    dg0 = fd.FiniteElement("DG", "triangle", 0)
    stress = fd.TensorFunctionSpace(mesh, dg0, symmetry=True)
    traction = fd.VectorFunctionSpace(mesh, dg0)
    mixed = velocity * stress * traction
    state = fd.Function(mixed)
    u, membrane, basal = fd.split(state)
    v, q, r = fd.TestFunctions(mixed)
    forcing = fd.as_vector((1.0, 1.0))
    if cubic_drag:
        # A uniform load has a uniform solution, on which every Krylov space
        # is one-dimensional: no solve could tell one Jacobian from another.
        x, y = fd.SpatialCoordinate(mesh)
        forcing = fd.as_vector((1.0 + 4.0 * fd.sin(6.0 * x), 1.0 - 3.0 * x * y))
    residual = (
        fd.inner(fd.grad(u), fd.grad(v))
        + (1.0 + cubic_drag * fd.inner(u, u)) * fd.inner(u, v)
        + fd.inner(
            membrane - fd.sym(fd.grad(u)),
            q - fd.sym(fd.grad(v)),
        )
        + fd.inner(basal - u, r - v)
        - fd.inner(forcing, v)
    ) * fd.dx
    # Match simulation.py's structural-zero blocks for retained-first SCPC.
    zero = fd.Constant(0.0)
    residual += fd.derivative(
        zero * membrane[0, 0] * basal[0] * fd.dx, state
    )
    return residual, state


def mixed_solver(mesh, mode, prefix, cubic_drag=0.0):
    """The solver simulation.py builds for ``mode``, on the toy residual."""
    residual, state = mixed_problem(mesh, cubic_drag)
    jacobian, pre_jacobian = None, None
    if linearization_state(mode) == "frozen":
        jacobian, pre_jacobian = frozen_linearization(residual, state)
    solver = fd.NonlinearVariationalSolver(
        fd.NonlinearVariationalProblem(residual, state, J=jacobian),
        solver_parameters=diagnostic_solver_parameters(mode),
        options_prefix=prefix,
        pre_jacobian_callback=pre_jacobian,
    )
    return solver, state


def test_mode(mesh, mode):
    os.environ["ISMIP7_DIAGNOSTIC_LINEAR_SOLVER"] = mode
    solver, state = mixed_solver(mesh, mode, f"ismip7_smoke_{mode}_")
    solver.solve()
    reason = solver.snes.getConvergedReason()
    if reason <= 0:
        raise RuntimeError(f"{mode} diverged with SNES reason {reason}")
    condensed = ""
    if mode.startswith("scpc_"):
        scpc = solver.snes.ksp.pc.getPythonContext()
        if scpc.condensed_solves < 1 or (
            scpc.condensed_iterations < scpc.condensed_solves
        ):
            raise RuntimeError(f"{mode} did not count its condensed solves")
        _, pmat = scpc.condensed_ksp.getOperators()
        near = pmat.getNearNullSpace()
        # An unset near-nullspace is a null handle; getVecs() on it segfaults.
        modes = near.getVecs() if near.handle else []
        expected = 3 if mode == "scpc_gamg" else 0
        if len(modes) != expected:
            raise RuntimeError(
                f"{mode} condensed operator has {len(modes)} near-nullspace "
                f"vectors, expected {expected}"
            )
        condensed = (
            f" condensed_solves={scpc.condensed_solves}"
            f" condensed_its={scpc.condensed_iterations}"
            f" near_nullspace={len(modes)}"
        )
    PETSc.Sys.Print(
        f"PASS {mode}: snes_its={solver.snes.getIterationNumber()} "
        f"linear_its={solver.snes.getLinearSolveIterations()}{condensed}"
    )
    return state


def test_frozen_linearization(mesh):
    r"""The line search's solve must see the Jacobian SCPC condensed.

    On a nonlinear residual the NLEQ-ERR line search solves for its simplified
    Newton step after evaluating the residual at a trial point. Condensed MUMPS
    is an exact inverse of the Jacobian at the Newton iterate, so with the
    linearization frozen there every outer solve takes one iteration; left
    live, the matrix-free operator has moved to the trial point and it cannot.
    Both must still reach the same root."""
    outcomes = {}
    for setting in ("1", "0"):
        os.environ["ISMIP7_FREEZE_LINEARIZATION"] = setting
        state_name = linearization_state("scpc_mumps")
        solver, state = mixed_solver(
            mesh, "scpc_mumps", f"ismip7_smoke_{state_name}_", cubic_drag=50.0
        )
        outer = []

        def monitor(ksp, its, rnorm, outer=outer):
            if its == 0:
                outer.append(0)
            else:
                outer[-1] = its

        solver.snes.ksp.setMonitor(monitor)
        solver.solve()
        newton = solver.snes.getIterationNumber()
        if newton < 3 or len(outer) <= newton:
            raise RuntimeError(
                f"{state_name}: not a nonlinear test ({newton} Newton "
                f"iterations, {len(outer)} linear solves)"
            )
        outcomes[state_name] = (state, newton, len(outer), max(outer), sum(outer))
    del os.environ["ISMIP7_FREEZE_LINEARIZATION"]

    frozen, live = outcomes["frozen"], outcomes["live"]
    if frozen[3] != 1:
        raise RuntimeError(
            f"frozen linearization: an exact condensed solve took {frozen[3]} "
            "outer iterations; the operator is not the one SCPC condensed"
        )
    if live[3] <= 1:
        raise RuntimeError("live linearization converged in one iteration: no test")
    difference = fd.errornorm(frozen[0], live[0]) / fd.norm(live[0])
    if difference > 1e-7:
        raise RuntimeError(f"frozen and live roots differ by {difference:.2e}")
    for name, (_, newton, solves, worst, total) in outcomes.items():
        PETSc.Sys.Print(
            f"PASS {name} linearization: snes_its={newton} linear_solves={solves} "
            f"outer_its_total={total} outer_its_max={worst}"
        )
    PETSc.Sys.Print(f"PASS frozen and live roots agree: rel diff {difference:.1e}")


def test_persistent_transport(mesh):
    space = fd.FunctionSpace(mesh, "DG", 0)
    thickness = fd.Function(space).assign(1.0)
    old = fd.Function(space)
    source = fd.Function(space)
    dt = fd.Constant(0.1)
    trial = fd.TrialFunction(space)
    test = fd.TestFunction(space)
    form = ((trial - old) / dt * test + trial * test - source * test) * fd.dx
    problem = fd.LinearVariationalProblem(
        fd.lhs(form), fd.rhs(form), thickness
    )
    solver = fd.LinearVariationalSolver(
        problem,
        solver_parameters=transport_solver_parameters(),
        options_prefix="ismip7_transport_smoke_",
    )
    expected = 1.0
    for value, timestep in ((1.0, 0.1), (2.0, 0.05)):
        old.assign(thickness)
        source.assign(value)
        dt.assign(timestep)
        solver.solve()
        ksp = solver.snes.getKSP()
        if ksp.getConvergedReason() <= 0:
            raise RuntimeError(
                "persistent transport solver diverged with KSP reason "
                f"{ksp.getConvergedReason()}"
            )
        expected = (expected + timestep * value) / (1.0 + timestep)
        target = fd.Function(space).assign(expected)
        if fd.errornorm(target, thickness) > 1e-12:
            raise RuntimeError("persistent transport solver used stale coefficients")
    PETSc.Sys.Print(
        "PASS persistent transport solver: source and dt coefficient updates; "
        f"last reason={ksp.getConvergedReason()} "
        f"iterations={ksp.getIterationNumber()} "
        f"residual={ksp.getResidualNorm():.3e}"
    )


def test_global_extrema(mesh):
    space = fd.FunctionSpace(mesh, "DG", 0)
    coordinates = fd.Function(
        fd.VectorFunctionSpace(mesh, space.ufl_element())
    ).interpolate(fd.SpatialCoordinate(mesh))
    x, y = fd.SpatialCoordinate(mesh)
    field = fd.Function(space).interpolate(x + 2.0 * y)
    lo, lo_xy = global_extreme_location(field, coordinates, mode="min")
    hi, hi_xy = global_extreme_location(field, coordinates, mode="max")
    if not (
        abs(lo - (lo_xy[0] + 2.0 * lo_xy[1])) < 1e-12
        and abs(hi - (hi_xy[0] + 2.0 * hi_xy[1])) < 1e-12
        and lo < hi
    ):
        raise RuntimeError("global extreme locations do not match field values")
    PETSc.Sys.Print(
        f"PASS global extrema: min={lo:.3f} at {lo_xy}, "
        f"max={hi:.3f} at {hi_xy}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "modes",
        nargs="*",
        default=["scpc_gamg", "scpc_mumps", "full_mumps"],
    )
    args = parser.parse_args()
    mesh = fd.UnitSquareMesh(2, 2)
    for mode in args.modes:
        test_mode(mesh, mode)
    test_frozen_linearization(fd.UnitSquareMesh(8, 8))
    test_persistent_transport(mesh)
    test_global_extrema(mesh)


if __name__ == "__main__":
    main()
