#!/usr/bin/env python
# coding: utf-8


import os
import firedrake
import icepack
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from icepack.constants import (
    ice_density as ρ_I,
    water_density as ρ_W,
    gravity as g,
)
from icepackaccs import extract_surface
from firedrake.petsc import PETSc
import datetime
import numpy as np

twostep = True

mesh_dir = "promice_meshes"
fn_template = mesh_dir + "/greenland_{:d}_{:d}.h5"

if mesh_dir == "promice_meshes":
    name = "promice"
    did = np.arange(2, 3153).tolist()
elif mesh_dir == "meshes":
    name = "detailed"
    did = np.arange(1, 2289).tolist()
else:
    name = "simple"
    did = np.arange(1, 64).tolist()

opts0 = {
    "dirichlet_ids": did,
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
        "snes_monitor": None,
        "snes_linesearch_monitor": None,
        # "ksp_monitor": None,
    },
}
opts1 = {
    "dirichlet_ids": did,
    # "side_wall_ids": [1, 3],
    "diagnostic_solver_type": "petsc",
    "diagnostic_solver_parameters": {
        "snes_type": "newtontr",
        "snes_tr_delta0": 1.0e5,
        "snes_tr_fallback_type": "dogleg",
        "snes_max_it": 5000,
        "snes_stol": 1.0e-8,
        "snes_rtol": 1.0e-8,
        "ksp_type": "bcgs",
        "ksp_max_it": 100000,
        "ksp_rtol": 1.0e-16,
        "ksp_atol": 1.0e-16,
        "pc_type": "bjacobi",
        "pc_hypre_type": "boomeramg",
        "pc_factor_mat_solver_type": "mumps",
        "pc_factor_shift_amount": 1.0e-10,
        "snes_monitor": None,
        "snes_linesearch_monitor": None,
        # "ksp_monitor": None,
    },
}

if firedrake.COMM_WORLD.size == 1:
    lcs = np.array([250 * 2 ** i for i in range(4, 8)])
elif firedrake.COMM_WORLD.size < 3:
    lcs = np.array([250 * 2 ** i for i in range(3, 8)])
elif firedrake.COMM_WORLD.size < 9:
    lcs = np.array([250 * 2 ** i for i in range(2, 6)])
elif firedrake.COMM_WORLD.size < 13:
    lcs = np.array([250 * 2 ** i for i in range(2, 6)])
else:
    lcs = np.array([250 * 2 ** i for i in range(0, 4)])

elapsed = np.zeros_like(lcs, dtype=float)
dampings = np.zeros_like(lcs, dtype=float)
delta0s = np.zeros_like(lcs, dtype=float)

timing_fn = "timings_{:s}_n{:d}.txt".format(name, firedrake.COMM_WORLD.size)
try_d0 = np.round(np.array([10 ** -(i / 2) * 1e7 for i in range(11)]), decimals=5)
try_damping = np.round(np.array([0.5 * 10 ** (-i / 2) for i in range(11)]), decimals=5)

with open(timing_fn, "w") as fout:
    fout.write("nproc, lc (m), time (s), LS damping, TR delta0\n")

for i, lc in enumerate(lcs):
    fn = fn_template.format(lc * 10, lc)
    with firedrake.CheckpointFile(fn, "r") as chk:
        mesh = chk.load_mesh(name="greenland")
        u_obs = chk.load_function(mesh, name="u")
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
    C_m = firedrake.assemble(C_var * firedrake.dx) / area
    C_mean = firedrake.Function(Q).interpolate(firedrake.Constant(C_m))
    C = C_var.copy(deepcopy=True)

    if firedrake.COMM_WORLD.size == 1:
        fig, ax = plt.subplots(2, 2, figsize=(14, 14))
        colors = firedrake.tripcolor(extract_surface(H), vmin=0, vmax=4000, axes=ax[0][0])
        plt.colorbar(colors, extend='max', label="Ice thickness [m]")
        colors = firedrake.tripcolor(extract_surface(u_obs_mag), cmap="turbo", norm=mcolors.LogNorm(vmin=5, vmax=2000), axes=ax[0][1])
        plt.colorbar(colors, extend='max', label='Velocity (m/yr)')
        colors = firedrake.tripcolor(extract_surface(C), vmin=0, vmax=0.5, cmap="plasma", axes=ax[1][0])
        plt.colorbar(colors, extend='max', label='C')
        colors = firedrake.tripcolor(extract_surface(C_var), vmin=0, vmax=0.5, cmap="plasma", axes=ax[1][1])
        plt.colorbar(colors, extend='max', label='C')
        fig.savefig(os.path.join("figs", os.path.split(os.path.splitext(fn)[0])[-1] + "_inputs.png"), dpi=300)

        for axeses in ax:
            for axes in axeses:
                axes.axis("equal")

    T = firedrake.Constant(260)
    A0 = icepack.rate_factor(T)

    def linear_pos_friction(**kwargs):
        p_W = ρ_W * g * firedrake.max_value(0, kwargs['thickness'] - kwargs['surface'])
        p_I = ρ_I * g * kwargs['thickness']
        ϕ = 1 - p_W / p_I
        return ϕ * kwargs['C'] ** 2.0 * firedrake.inner(kwargs['velocity'], kwargs['velocity'])

    model = icepack.models.HybridModel(friction=linear_pos_friction)

    PETSc.Sys.Print("Beginning initial velocity solve {:d}".format(lc))

    u0 = firedrake.Function(V).interpolate(firedrake.Constant(1.0e-1) * u_obs)
    if twostep:
        for damping in try_damping:
            opts0["diagnostic_solver_parameters"]["snes_linesearch_damping"] = damping
            solver0 = icepack.solvers.FlowSolver(model, **opts0)
            PETSc.Sys.Print('Trying LS solve with damping = {:e}'.format(damping))
            try:
                start_time = datetime.datetime.now()
                u0 = solver0.diagnostic_solve(
                    velocity=u0,
                    fluidity=A0,
                    C=C,
                    surface=h,
                    thickness=H
                )
                end_time = datetime.datetime.now()
                elapsed[i] = (end_time - start_time).seconds + (end_time - start_time).microseconds * 1e-6
                dampings[i] = damping
                break
            except firedrake.exceptions.ConvergenceError:
                pass
        else:
            dampings[i] = np.nan
            PETSc.Sys.Print('LS Failed')
    PETSc.Sys.Print('Starting TR solve on {:d} processes'.format(firedrake.COMM_WORLD.size))
    for d0 in try_d0:
        opts1["diagnostic_solver_parameters"]["snes_tr_delta0"] = d0
        solver1 = icepack.solvers.FlowSolver(model, **opts1)
        PETSc.Sys.Print('Trying TR solve with d0 = {:e}'.format(d0))
        try:
            start_time = datetime.datetime.now()
            u0 = solver1.diagnostic_solve(
                velocity=u0,
                fluidity=A0,
                C=C,
                surface=h,
                thickness=H
            )
            end_time = datetime.datetime.now()
            elapsed[i] += (end_time - start_time).seconds + (end_time - start_time).microseconds * 1e-6
            delta0s[i] = d0
            break
        except firedrake.exceptions.ConvergenceError:
            pass
    else:
        PETSc.Sys.Print('No solution found!')
        delta0s[i] = np.nan

    PETSc.Sys.Print('Solve on {:d} processes took {:d} s'.format(firedrake.COMM_WORLD.size, int(elapsed[i])))

    with open(timing_fn, "a") as fout:
        fout.write("{:d}, {:d}, {:f}, {:e}, {:e}\n".format(firedrake.COMM_WORLD.size, lc, elapsed[i], dampings[i], delta0s[i]))

    if firedrake.COMM_WORLD.size == 1:
        fig, ax = plt.subplots(1, 2, figsize=(14, 7))
        colors = firedrake.tripcolor(extract_surface(u_obs_mag), vmin=0, vmax=1000, axes=ax[0])
        plt.colorbar(colors, extend='max', label="Observed flow speed [m/yr]")
        colors = firedrake.tripcolor(extract_surface(firedrake.sqrt(u0[0] ** 2.0 + u0[1] ** 2.0)), vmin=0, vmax=1000, axes=ax[1])
        plt.colorbar(colors, extend='max', label="Modeled flow speed [m/yr]")
        ax[0].axis("equal")
        ax[1].axis("equal")
        fig.savefig(os.path.join("figs", os.path.split(os.path.splitext(fn)[0])[-1] + "_veltest.png"), dpi=300)
