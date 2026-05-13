#! /usr/bin/env python3
# vim:fenc=utf-8
#
# Copyright © 2026 David Lilien <dlilien@iu.edu>
#
"""
Run inversion tests on initialized Greenland meshes.

Usage:
    inversion_tests.py --mesh-kind <mesh_kind>

    inversion_tests.py --help
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_MESH_DIR = Path(__file__).resolve().parent.parent / "meshes"
DEFAULT_LOG_DIR = Path("test_logs")
MESH_KINDS = ("detailed", "promice", "simple", "buffered")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run inversion tests.")
    parser.add_argument(
        "--mesh-dir",
        type=Path,
        default=DEFAULT_MESH_DIR,
        help="Directory containing initialized h5 files.",
    )
    parser.add_argument(
        "--mesh-kind",
        choices=MESH_KINDS,
        default="detailed",
        help="Mesh type to use: detailed, promice, simple, or buffered.",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help="Directory for timing logs and cached inversion results.",
    )
    parser.add_argument(
        "--lambda-exponent",
        type=int,
        default=6,
        help="Use 10 ** (2 * lambda_exponent) for the regularization length.",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=3,
        help="Maximum optimization iterations.",
    )
    return parser.parse_args()


def dirichlet_ids(mesh_kind: str) -> list[int]:
    """Return Dirichlet boundary IDs for a mesh kind."""
    if mesh_kind == "promice":
        return np.arange(2, 3153).tolist()
    if mesh_kind == "detailed":
        return np.arange(1, 2289).tolist()
    return np.arange(1, 64).tolist()


def solver_options(dids: list[int]) -> tuple[dict, dict]:
    """Build line-search and trust-region solver option dictionaries."""
    opts0 = {
        "dirichlet_ids": dids,
        "diagnostic_solver_type": "petsc",
        "diagnostic_solver_parameters": {
            "snes_type": "newtonls",
            "snes_line_search_type": "bt",
            "snes_linesearch_order": 2,
            "snes_linesearch_max_it": 2500,
            "snes_linesearch_damping": 0.05,
            "snes_max_it": 5000,
            "snes_stol": 1.0e-6,
            "snes_rtol": 1.0e-5,
            "ksp_type": "bcgs",
            "ksp_max_it": 2500,
            "ksp_rtol": 1.0e-8,
            "ksp_atol": 1.0e-4,
            "ksp_converged_maxits": True,
            "pc_type": "bjacobi",
            "pc_factor_mat_solver_type": "mumps",
            "pc_factor_shift_amount": 1.0e-10,
        },
    }
    opts1 = {
        "dirichlet_ids": dids,
        "diagnostic_solver_type": "petsc",
        "diagnostic_solver_parameters": {
            "snes_type": "newtontr",
            "snes_tr_delta0": 1.0e5,
            "snes_tr_fallback_type": "dogleg",
            "snes_max_it": 5000,
            "snes_stol": 1.0e-8,
            "snes_rtol": 1.0e-8,
            "snes_atol": 1.0e-3,
            "ksp_type": "bcgs",
            "ksp_max_it": 100000,
            "ksp_rtol": 1.0e-16,
            "ksp_atol": 1.0e-16,
            "pc_type": "bjacobi",
            "pc_hypre_type": "boomeramg",
            "pc_factor_mat_solver_type": "mumps",
            "pc_factor_shift_amount": 1.0e-10,
        },
    }
    return opts0, opts1


def all_resolution_levels() -> np.ndarray:
    """Return all mesh resolution levels."""
    return np.array([250 * 2**i for i in range(0, 8)])


def selected_resolution_levels(comm_size: int) -> np.ndarray:
    """Choose inversion test resolutions based on MPI size."""
    all_lcs = all_resolution_levels()
    if comm_size == 1:
        return all_lcs[8:3:-1]
    if comm_size < 3:
        return all_lcs[6:2:-1]
    if comm_size < 9:
        return all_lcs[5:1:-1]
    if comm_size < 13:
        return all_lcs[2:6]
    return all_lcs[:4]


def checkpoint_path(mesh_dir: Path, mesh_kind: str, lc: int) -> Path:
    """Build the checkpoint path matching mesh_greenland.py names."""
    return mesh_dir / f"greenland_{mesh_kind}_{10 * lc}_{lc}.h5"


def timing_input_path(log_dir: Path, mesh_kind: str, comm_size: int) -> Path:
    """Find diagnostic timing input from test_logs or the current directory."""
    timing_fn = f"timings_{mesh_kind}_n{comm_size:d}.txt"
    for path in (log_dir / timing_fn, Path(timing_fn)):
        if path.exists():
            return path
    raise FileNotFoundError(
        f"Could not find diagnostic timing file {timing_fn} in {log_dir} or cwd"
    )


def set_solver_parameters_from_timings(
    opts0: dict, opts1: dict, timing_in, lc: int
) -> None:
    """Set solver parameters from diagnostic timing results."""
    rows = timing_in[timing_in["lc (m)"] == lc]
    if len(rows) == 0:
        raise ValueError(f"No diagnostic timing row for lc={lc}")
    opts0["diagnostic_solver_parameters"]["snes_linesearch_damping"] = float(
        rows["LS damping"].values[0]
    )
    opts1["diagnostic_solver_parameters"]["snes_tr_delta0"] = float(
        rows["TR delta0"].values[0]
    )


def run_inversion(
    mesh_dir: Path,
    mesh_kind: str,
    log_dir: Path,
    lambda_exponent: int,
    max_iterations: int,
) -> None:
    """Run inversion tests for one mesh kind."""
    import firedrake
    import icepack
    from firedrake.petsc import PETSc
    from icepack.constants import gravity as g
    from icepack.constants import ice_density as ρ_I
    from icepack.constants import water_density as ρ_W
    from icepack.statistics import MaximumProbabilityEstimator
    from icepack.statistics import StatisticsProblem

    log_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = log_dir / "cached_results"
    cache_dir.mkdir(parents=True, exist_ok=True)

    opts0, opts1 = solver_options(dirichlet_ids(mesh_kind))
    lcs = selected_resolution_levels(firedrake.COMM_WORLD.size)
    all_lcs = all_resolution_levels()

    timing_in = pd.read_csv(
        timing_input_path(log_dir, mesh_kind, firedrake.COMM_WORLD.size),
        delimiter=", ",
        engine="python",
    )
    timing_out_fn = (
        log_dir / f"inv_timings_{mesh_kind}_n{firedrake.COMM_WORLD.size:d}.txt"
    )
    out_template = f"inv_{mesh_kind}_{{:d}}.h5"
    LCap = 10 ** (2 * lambda_exponent)

    with timing_out_fn.open("w") as fout:
        fout.write("nproc, lc (m), lambda, max iterations\n")

    for lc in lcs:
        lc = int(lc)
        set_solver_parameters_from_timings(opts0, opts1, timing_in, lc)

        checkpoint_fn = checkpoint_path(mesh_dir, mesh_kind, lc)
        if not checkpoint_fn.exists():
            print(f"Skipping missing checkpoint: {checkpoint_fn}")
            continue

        with firedrake.CheckpointFile(str(checkpoint_fn), "r") as chk:
            mesh = chk.load_mesh(name="greenland")
            u_obs = chk.load_function(mesh, name="u")
            sigma_u = chk.load_function(mesh, name="sigma_u")
            H_in = chk.load_function(mesh, name="H")
            b_in = chk.load_function(mesh, name="b")

        area = firedrake.Constant(
            firedrake.assemble(firedrake.Constant(1.0) * firedrake.ds_t(mesh))
        )
        Q = firedrake.FunctionSpace(mesh, "CG", 2, vfamily="R", vdegree=0)
        V = firedrake.VectorFunctionSpace(mesh, "CG", 1, dim=2, vfamily="GL", vdegree=4)

        h0 = firedrake.Function(Q).interpolate(b_in + H_in)
        H = firedrake.Function(Q).interpolate(H_in)
        h = h0.copy(deepcopy=True)

        # Smooth the surface elevation to reduce bumps etc.
        α = firedrake.Constant(2e3)
        J = (
            0.5
            * (
                (h - h0) ** 2
                + α**2 * firedrake.inner(firedrake.grad(h), firedrake.grad(h))
            )
            * firedrake.dx
        )
        F = firedrake.derivative(J, h)
        firedrake.solve(F == 0, h)

        τ_D = firedrake.Function(Q).interpolate(
            ρ_I
            * g
            * H
            * firedrake.sqrt(firedrake.grad(H)[0] ** 2.0 + firedrake.grad(H)[1] ** 2.0)
        )
        u_obs_mag = firedrake.Function(Q).interpolate(
            firedrake.max_value(
                firedrake.sqrt(u_obs[0] ** 2 + u_obs[1] ** 2),
                firedrake.Constant(1.0e-3),
            )
        )
        C_var = firedrake.Function(Q).interpolate(firedrake.sqrt(τ_D / u_obs_mag / 2.0))
        C = C_var.copy(deepcopy=True)

        T = firedrake.Constant(260)
        A0 = icepack.rate_factor(T)

        def linear_pos_friction(**kwargs):
            """Evaluate flotation-adjusted linear sliding friction."""
            p_W = (
                ρ_W
                * g
                * firedrake.max_value(0, kwargs["thickness"] - kwargs["surface"])
            )
            p_I = ρ_I * g * kwargs["thickness"]
            ϕ = 1 - p_W / p_I
            return (
                ϕ
                * kwargs["C"] ** 2.0
                * firedrake.inner(kwargs["velocity"], kwargs["velocity"])
            )

        model = icepack.models.HybridModel(friction=linear_pos_friction)
        solver0 = icepack.solvers.FlowSolver(model, **opts0)
        solver1 = icepack.solvers.FlowSolver(model, **opts1)

        PETSc.Sys.Print(f"Beginning initial velocity solve {lc:d}")
        u0 = firedrake.Function(V).interpolate(firedrake.Constant(1.0e-1) * u_obs)
        u0 = solver0.diagnostic_solve(
            velocity=u0, fluidity=A0, C=C, surface=h, thickness=H
        )

        targ_index = np.where(all_lcs == lc)[0][0]
        if targ_index < (len(all_lcs) - 1):
            lowres_cache_fn = cache_dir / out_template.format(
                int(all_lcs[targ_index + 1])
            )
            if not lowres_cache_fn.exists():
                raise FileNotFoundError(f"Run lowres first: {lowres_cache_fn}")
            with firedrake.CheckpointFile(str(lowres_cache_fn), "r") as chk:
                lowres_mesh = chk.load_mesh("greenland")
                C_old = chk.load_function(lowres_mesh, "C")
                C = firedrake.Function(Q).interpolate(C_old)
        else:
            C = C_var.copy(deepcopy=True)
        u = u0.copy(deepcopy=True)

        def simulation(C):
            """Run the diagnostic model for the supplied friction field."""
            return solver1.diagnostic_solve(
                velocity=u, fluidity=A0, C=C, surface=h, thickness=H
            )

        def loss_functional(u):
            """Evaluate the normalized velocity misfit functional."""
            δu = u - u_obs
            return (
                0.5
                / area
                * ((δu[0] / sigma_u[0]) ** 2 + (δu[1] / sigma_u[1]) ** 2)
                * firedrake.ds_t(mesh)
            )

        def total_misfit(u):
            """Evaluate the mean speed misfit functional."""
            δu = u - u_obs
            return (
                0.5
                / area
                * firedrake.sqrt((δu[0]) ** 2 + (δu[1]) ** 2)
                * firedrake.ds_t(mesh)
            )

        def regularization(C):
            """Evaluate the friction-field smoothing regularization."""
            L = firedrake.Constant(LCap)
            return (
                0.5
                / area
                * L**2
                * firedrake.inner(firedrake.grad(C), firedrake.grad(C))
                * firedrake.ds_b(mesh)
            )

        C = firedrake.Function(Q).interpolate(C)
        u = simulation(C)
        problem = StatisticsProblem(
            simulation=simulation,
            loss_functional=loss_functional,
            regularization=regularization,
            controls=C,
        )
        estimator = MaximumProbabilityEstimator(
            problem,
            gradient_tolerance=1e-4,
            step_tolerance=1e-4,
            max_iterations=max_iterations,
        )
        PETSc.Sys.Print("Optimizing", out_template.format(lc))

        C = estimator.solve()
        u_opt = simulation(C)
        state = estimator._solver.getAlgorithmState()

        cache_fn = cache_dir / out_template.format(lc)
        with firedrake.CheckpointFile(str(cache_fn), "w") as chk:
            chk.create_group("metadata")
            chk.set_attr("metadata", "lambda", LCap)
            chk.set_attr("metadata", "max_iterations", max_iterations)
            chk.set_attr("metadata", "gnorm", state.gnorm)
            chk.set_attr("metadata", "cnorm", state.cnorm)
            chk.set_attr("metadata", "snorm", state.snorm)
            chk.set_attr("metadata", "loss", firedrake.assemble(loss_functional(u)))
            chk.set_attr(
                "metadata", "regularization", firedrake.assemble(regularization(C))
            )
            chk.set_attr(
                "metadata",
                "cost",
                firedrake.assemble(loss_functional(u))
                + firedrake.assemble(regularization(C)),
            )
            chk.set_attr(
                "metadata", "average_misfit", firedrake.assemble(total_misfit(u))
            )
            chk.save_mesh(mesh)
            chk.save_function(C, name="C")
            chk.save_function(u_opt, name="u_opt")
            chk.save_function(u_obs, name="u_obs")

        with timing_out_fn.open("a") as fout:
            fout.write(
                f"{firedrake.COMM_WORLD.size:d}, {lc:d}, {LCap:e}, {max_iterations:d}\n"
            )


def main() -> None:
    """Run inversion tests from command-line arguments."""
    args = parse_args()
    sys.argv = [sys.argv[0]]
    run_inversion(
        args.mesh_dir.expanduser(),
        args.mesh_kind,
        args.log_dir.expanduser(),
        args.lambda_exponent,
        args.max_iterations,
    )


if __name__ == "__main__":
    main()
