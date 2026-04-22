#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import firedrake
import icepack
from icepack.constants import (
    ice_density as ρ_I,
    water_density as ρ_W,
    gravity as g,
)
import numpy as np
from icepack.statistics import (
    StatisticsProblem,
    MaximumProbabilityEstimator,
)

from firedrake.petsc import PETSc
import pandas as pd

Lexp = 6

mesh_dir = "meshes"
fn_template = mesh_dir + "/greenland_{:d}_{:d}.h5"

if mesh_dir == "promice_meshes":
    name = "promice"
    dirichlet_ids = np.arange(2, 3153).tolist()
if mesh_dir == "meshes":
    name = "detailed"
    dirichlet_ids = np.arange(1, 2289).tolist()
else:
    name = "simple"
    dirichlet_ids = np.arange(1, 64).tolist()

opts0 = {
    "dirichlet_ids": dirichlet_ids,
    # "side_wall_ids": [1, 3],
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
        # "pc_hypre_type": "boomeramg",
        "pc_factor_mat_solver_type": "mumps",
        "pc_factor_shift_amount": 1.0e-10,
    },
}
opts1 = {
    "dirichlet_ids": dirichlet_ids,
    # "side_wall_ids": [1, 3],
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

all_lcs = np.array([250 * 2 ** i for i in range(0, 8)])
if firedrake.COMM_WORLD.size == 1:
    lcs = all_lcs[8:3:-1]
elif firedrake.COMM_WORLD.size < 3:
    lcs = all_lcs[6:2:-1]
elif firedrake.COMM_WORLD.size < 9:
    lcs = all_lcs[5:1:-1]
elif firedrake.COMM_WORLD.size < 13:
    lcs = all_lcs[2:6]
else:
    lcs = all_lcs[:4]


elapsed = np.zeros_like(lcs, dtype=float)

timing_in_fn = "timings_{:s}_n{:d}.txt".format(name, firedrake.COMM_WORLD.size)
timing_out_fn = "inv_timings_{:s}_n{:d}.txt".format(name, firedrake.COMM_WORLD.size)
out_template = "inv_{:s}_{:d}.h5"

timing_in = pd.read_csv(timing_in_fn, delimiter=", ", engine="python")


for lc in lcs:
    opts0["diagnostic_solver_parameters"]["snes_linesearch_damping"] = float(timing_in[timing_in["lc (m)"] == lc]["LS damping"].values[0])
    opts1["diagnostic_solver_parameters"]["snes_tr_delta0"] = float(timing_in[timing_in["lc (m)"] == lc]["TR delta0"].values[0])

    fn = fn_template.format(lc * 10, lc)
    with firedrake.CheckpointFile(fn, "r") as chk:
        mesh = chk.load_mesh(name="greenland")
        u_obs = chk.load_function(mesh, name="u")
        sigma_u = chk.load_function(mesh, name="sigma_u")
        smb = chk.load_function(mesh, name="smb")
        H_in = chk.load_function(mesh, name="H")
        b_in = chk.load_function(mesh, name="b")


    area = firedrake.Constant(firedrake.assemble(firedrake.Constant(1.0) * firedrake.ds_t(mesh)))
    Q = firedrake.FunctionSpace(
        mesh, "CG", 2, vfamily="R", vdegree=0
    )
    V = firedrake.VectorFunctionSpace(
        mesh, "CG", 1, dim=2, vfamily="GL", vdegree=4
    )
    x, y, ζ = firedrake.SpatialCoordinate(mesh)

    h0 = firedrake.Function(Q).interpolate(b_in + H_in)
    H0 = firedrake.Function(Q).interpolate(H_in)
    b = firedrake.Function(Q).interpolate(b_in)

    h = h0.copy(deepcopy=True)
    H = H0.copy(deepcopy=True)

    # Smooth the surface elevation to reduce bumps etc.
    α = firedrake.Constant(2e3)
    J = 0.5 * ((h - h0) ** 2 + α ** 2 * firedrake.inner(firedrake.grad(h), firedrake.grad(h))) * firedrake.dx
    F = firedrake.derivative(J, h)
    firedrake.solve(F == 0, h)

    τ_D = firedrake.Function(Q).interpolate(ρ_I * g * H * firedrake.sqrt(firedrake.grad(H)[0] ** 2.0 + firedrake.grad(H)[1] ** 2.0))

    u_obs_mag = firedrake.Function(Q).interpolate(firedrake.max_value(firedrake.sqrt(u_obs[0]**2 + u_obs[1]**2), firedrake.Constant(1.0e-3)))

    C_var = firedrake.Function(Q).interpolate(firedrake.sqrt(τ_D / u_obs_mag / 2.0))
    C = C_var.copy(deepcopy=True)

    T = firedrake.Constant(260)
    A0 = icepack.rate_factor(T)

    def linear_pos_friction(**kwargs):
        p_W = ρ_W * g * firedrake.max_value(0, kwargs['thickness'] - kwargs['surface'])
        p_I = ρ_I * g * kwargs['thickness']
        ϕ = 1 - p_W / p_I
        return ϕ * kwargs['C'] ** 2.0 * firedrake.inner(kwargs['velocity'], kwargs['velocity'])

    model = icepack.models.HybridModel(friction=linear_pos_friction)

    solver0 = icepack.solvers.FlowSolver(model, **opts0)
    solver1 = icepack.solvers.FlowSolver(model, **opts1)

    PETSc.Sys.Print("Beginning initial velocity solve {:d}".format(lc))

    u0 = firedrake.Function(V).interpolate(firedrake.Constant(1.0e-1) * u_obs)
    σ = firedrake.Function(Q).interpolate(firedrake.sqrt(sigma_u[0]**2 + sigma_u[1]**2))

    # Ballpark correct solution, but not too precise so that we do not have convergence issues due to no change
    u0 = solver0.diagnostic_solve(
        velocity=u0,
        fluidity=A0,
        C=C,
        surface=h,
        thickness=H
    )


    LCap = 10 ** (2 * Lexp)
    targ_index = np.where(all_lcs == lc)[0][0]
    if targ_index < (len(all_lcs) - 1):
        lowres_cache_fn = os.path.join("cached_results", out_template.format(name, all_lcs[targ_index + 1]))
        if not os.path.exists(lowres_cache_fn):
            raise FileNotFoundError("Run lowres first")
        with firedrake.CheckpointFile(lowres_cache_fn, "r") as chk:
            lowres_mesh = chk.load_mesh("greenland")
            C_old = chk.load_function(lowres_mesh, "C")
            C = firedrake.Function(Q).interpolate(C_old)
    else:
        C = C_var.copy(deepcopy=True)
    u = u0.copy(deepcopy=True)

    def simulation(C):
        return solver1.diagnostic_solve(
            velocity=u,
            fluidity=A0,
            C=C,
            surface=h,
            thickness=H
        )

    def loss_functional(u):
        δu = u - u_obs
        return 0.5 / area * ((δu[0] / sigma_u[0])**2 + (δu[1] / sigma_u[1])**2) * firedrake.ds_t(mesh)

    def total_misfit(u):
        δu = u - u_obs
        return 0.5 / area * firedrake.sqrt((δu[0])**2 + (δu[1])**2) * firedrake.ds_t(mesh)

    def regularization(C):
        L = firedrake.Constant(LCap)
        return 0.5 / area * (L)**2 * firedrake.inner(firedrake.grad(C), firedrake.grad(C)) * firedrake.ds_b(mesh)

    C = firedrake.Function(Q).interpolate(C)
    u = simulation(C)

    problem = StatisticsProblem(
        simulation=simulation,
        loss_functional=loss_functional,
        regularization=regularization,
        controls=C,
    )

    n_iter = 3
    estimator = MaximumProbabilityEstimator(
        problem,
        gradient_tolerance=1e-4,
        step_tolerance=1e-4,
        max_iterations=n_iter,
    )
    PETSc.Sys.Print("Optimizing", out_template.format(name, lc))

    C = estimator.solve()

    u_opt = simulation(C)

    state = estimator._solver.getAlgorithmState()
    cache_fn = os.path.join("cached_results", out_template.format(name, lc))
    with firedrake.CheckpointFile(cache_fn, "w") as chk:
        chk.create_group("metadata")
        chk.set_attr("metadata", "lambda", LCap)
        # chk.set_attr("metadata", "iterations", n_iter_elapsed)
        chk.set_attr("metadata", "max_iterations", n_iter)
        chk.set_attr("metadata", "gnorm", state.gnorm)
        chk.set_attr("metadata", "cnorm", state.cnorm)
        chk.set_attr("metadata", "snorm", state.snorm)
        chk.set_attr("metadata", "loss", firedrake.assemble(loss_functional(u)))
        chk.set_attr("metadata", "regularization", firedrake.assemble(regularization(C)))
        chk.set_attr("metadata", "cost", firedrake.assemble(loss_functional(u)) + firedrake.assemble(regularization(C)))
        chk.set_attr("metadata", "average_misfit", firedrake.assemble(total_misfit(u)))

        chk.save_mesh(mesh)
        chk.save_function(C, name="C")
        chk.save_function(u_opt, name="u_opt")
        chk.save_function(u_obs, name="u_obs")
