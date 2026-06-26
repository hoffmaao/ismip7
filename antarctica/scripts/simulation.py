#!/usr/bin/env python3
r"""Shared simulation engine for ISMIP7 Antarctic experiments."""

import numpy as np
import os, sys, glob, json
from time import perf_counter

import firedrake as fd
from firedrake import (
    Constant,
    Function,
    max_value,
    sqrt,
    inner,
    grad,
    derivative,
    dx,
    dS,
    ds,
    split,
    assemble,
    Mesh,
    FunctionSpace,
    VectorFunctionSpace,
    TensorFunctionSpace,
    FiniteElement,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    exp,
)
from firedrake.petsc import PETSc

import rasterio, icepack
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
    glen_flow_law as n_glen_val,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
RESULTS_DIR = os.path.join(_ROOT, "results")

lc = int(os.environ.get("ISMIP7_LC", "2500"))
lc_coarse = int(os.environ.get("ISMIP7_LC_COARSE", "64000"))

# SI ice density (kg/m^3) for VAF/mass diagnostics only. The dynamics use the
# icepack rho_I (MPa-m-yr units); thickness fields are lengths (m) in both, so
# the Gt/SLE conversion just needs the SI density.
RHO_I_SI = 917.0
GT_PER_MM_SLE = 362.5  # Gt of water per mm global-mean sea-level rise


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def setup_model(restart_from=None):
    r"""Load mesh, data, inversion fields, and build diagnostic solver."""
    os.makedirs(RESULTS_DIR, exist_ok=True)

    mesh_fn = os.environ.get(
        "ISMIP7_MESH", os.path.join(MESH_DIR, f"antarctica_{lc_coarse}_{lc}.msh")
    )
    PETSc.Sys.Print(f"Loading mesh: {mesh_fn}")
    mesh = Mesh(mesh_fn)
    PETSc.Sys.Print(f"  {mesh.num_vertices()} vertices, {mesh.num_cells()} cells")

    bndids_fn = os.environ.get(
        "ISMIP7_BNDIDS", os.path.join(MESH_DIR, "boundary_ids.json")
    )
    with open(bndids_fn) as f:
        bnd_ids = json.load(f)
    calving_ids = tuple(bnd_ids["calving"])
    use_calving_terminus = os.environ.get("ISMIP7_NO_CALVING_TERMINUS") is None

    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Sigma = TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T

    PETSc.Sys.Print("Loading data...")
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    b = icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:bed"), Q)
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "10.0"))
    H = Function(Q, name="thickness").interpolate(
        max_value(
            icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q),
            Constant(h_clamp),
        )
    )
    rho_ratio = Constant(917.0 / 1024.0)
    s = Function(Q, name="surface").interpolate(
        max_value(b + H, (Constant(1.0) - rho_ratio) * H)
    )

    vel_fn = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:VX"), rasterio.open(f"netcdf:{vel_fn}:VY")),
        V,
        fillvalue=0.0,
    )

    inv_fn = os.path.join(MESH_DIR, f"inversion_icepack2_{lc}.h5")
    PETSc.Sys.Print(f"Loading MAP: {inv_fn}")
    with fd.CheckpointFile(inv_fn, "r") as chk:
        chk_mesh = chk.load_mesh()
        theta = chk.load_function(chk_mesh, name="log_friction")
        phi = chk.load_function(chk_mesh, name="log_fluidity")
    theta_f = Function(Q, name="theta")
    theta_f.dat.data[:] = theta.dat.data_ro
    phi_f = Function(Q, name="phi")
    phi_f.dat.data[:] = phi.dat.data_ro

    A0 = Constant(icepack.rate_factor(Constant(260.0)))
    n = Constant(n_glen_val)
    tau_c = Constant(0.1)
    u_c = Constant(float(Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    ).dat.data_ro.mean()))

    phi_eff = Function(Q).interpolate(
        max_value(
            Constant(1.0)
            - rho_W * g * max_value(Constant(0.0), -b) / (rho_I * g * max_value(H, Constant(1.0))),
            Constant(0.01),
        )
    )
    K_base = u_c / (phi_eff * tau_c) ** n
    A_map = A0 * exp(phi_f)
    K_map = K_base * exp(-n * theta_f)

    sparams = {
        "snes_type": "newtonls",
        "snes_max_it": 200,
        "snes_linesearch_type": "nleqerr",
        "snes_divergence_tolerance": -1,
        "snes_stol": 0.0,
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 400,
        "mat_mumps_icntl_24": 1,
        "mat_mumps_cntl_3": 1e-12,
    }
    fc_params = {"quadrature_degree": 4}

    z = Function(Z)
    z.sub(0).interpolate(Constant(0.1) * u_obs)

    h = H.copy(deepcopy=True)
    h.rename("thickness")

    if restart_from is not None:
        PETSc.Sys.Print(f"Restarting from: {restart_from}")
        with fd.CheckpointFile(restart_from, "r") as chk:
            rst_mesh = chk.load_mesh()
            h_rst = chk.load_function(rst_mesh, name="thickness")
            u_rst = chk.load_function(rst_mesh, name="velocity")
        h.dat.data[:] = h_rst.dat.data_ro
        z.sub(0).dat.data[:] = u_rst.dat.data_ro
        s.interpolate(max_value(b + h, (Constant(1.0) - rho_ratio) * h))
        phi_eff.interpolate(
            max_value(
                Constant(1.0)
                - rho_W * g * max_value(Constant(0.0), -b)
                / (rho_I * g * max_value(h, Constant(1.0))),
                Constant(0.01),
            )
        )

    u_s, M_s, tau_s = split(z)
    fields = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": h,
        "surface": s,
    }
    rheology = {
        "flow_law_exponent": n,
        "flow_law_coefficient": A_map,
        "sliding_exponent": n,
        "sliding_coefficient": K_map,
    }

    L = (
        model.minimization.viscous_power(**fields, **rheology)
        + model.minimization.friction_power(**fields, **rheology)
        + model.minimization.momentum_balance(**fields)
    )
    if use_calving_terminus:
        L += model.minimization.calving_terminus(**fields, outflow_ids=calving_ids)

    prob = NonlinearVariationalProblem(
        derivative(L, z), z, form_compiler_parameters=fc_params
    )
    slvr = NonlinearVariationalSolver(prob, solver_parameters=sparams)

    PETSc.Sys.Print("Initial diagnostic solve...")
    for exponent in np.linspace(1.0, n_glen_val, 5):
        n.assign(exponent)
        slvr.solve()
    PETSc.Sys.Print("  Done")

    accum = Function(Q, name="accumulation").assign(0.0)

    return {
        "mesh": mesh,
        "Q": Q,
        "V": V,
        "Z": Z,
        "z": z,
        "h": h,
        "s": s,
        "b": b,
        "slvr": slvr,
        "n": n,
        "accum": accum,
        "phi_eff": phi_eff,
        "rho_ratio": rho_ratio,
        "h_clamp": h_clamp,
        "calving_ids": calving_ids,
        "u_obs": u_obs,
    }


def run_simulation(
    ctx,
    experiment_name,
    t_start,
    t_end,
    dt=1.0,
    output_interval=10,
    checkpoint_interval=100,
    forcing_callback=None,
):
    r"""Run the split diagnostic-prognostic time-stepping loop."""
    mesh = ctx["mesh"]
    Q = ctx["Q"]
    z = ctx["z"]
    h = ctx["h"]
    s = ctx["s"]
    b = ctx["b"]
    slvr = ctx["slvr"]
    n = ctx["n"]
    accum = ctx["accum"]
    phi_eff = ctx["phi_eff"]
    rho_ratio = ctx["rho_ratio"]
    h_clamp = ctx["h_clamp"]

    nsteps = int((t_end - t_start) / dt)
    dt_c = Constant(dt)
    PETSc.Sys.Print(
        f"\nTime-stepping: {t_start}->{t_end}, dt={dt}yr, {nsteps} steps"
    )

    Q_dg = FunctionSpace(mesh, "DG", 0)
    h_dg = Function(Q_dg, name="h_dg")
    h_dg_old = Function(Q_dg)
    phi_dg = fd.TestFunction(Q_dg)
    h_dg_trial = fd.TrialFunction(Q_dg)

    n_facet = fd.FacetNormal(mesh)

    s_float = Function(Q).interpolate(
        b + (rho_W / rho_I) * max_value(-b, Constant(0.0))
    )

    results = []

    for k in range(1, nsteps + 1):
        t_step_start = perf_counter()
        t_yr = t_start + k * dt

        if forcing_callback is not None:
            forcing_callback(ctx, t_yr)

        try:
            slvr.solve()
        except fd.ConvergenceError:
            PETSc.Sys.Print(f"  Step {k}: re-doing continuation...")
            for exponent in np.linspace(1.0, n_glen_val, 5):
                n.assign(exponent)
                slvr.solve()

        u_vel = z.subfunctions[0]

        # DG0 upwind implicit Euler for thickness
        h_dg.project(h)
        h_dg_old.assign(h_dg)

        un = fd.dot(u_vel, n_facet)
        un_plus = (un + abs(un)) / 2

        F_prog = (
            (h_dg_trial - h_dg_old) / dt_c * phi_dg * dx
            - h_dg_trial * fd.div(u_vel * phi_dg) * dx
            + (un_plus("+") * h_dg_trial("+") - un_plus("-") * h_dg_trial("-"))
            * fd.jump(phi_dg)
            * dS
            + un_plus * h_dg_trial * phi_dg * ds
            - accum * phi_dg * dx
        )
        fd.solve(
            fd.lhs(F_prog) == fd.rhs(F_prog),
            h_dg,
            solver_parameters={
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",
            },
        )

        h_dg.interpolate(max_value(h_dg, Constant(h_clamp)))
        h.project(h_dg)
        h.interpolate(max_value(h, Constant(h_clamp)))

        s.interpolate(max_value(b + h, (Constant(1.0) - rho_ratio) * h))

        phi_eff.interpolate(
            max_value(
                Constant(1.0)
                - rho_W * g * max_value(Constant(0.0), -b)
                / (rho_I * g * max_value(h, Constant(1.0))),
                Constant(0.01),
            )
        )

        t_elapsed = perf_counter() - t_step_start

        haf = Function(Q).interpolate(max_value(s - s_float, Constant(0.0)))
        vaf = float(assemble(haf * dx)) * RHO_I_SI / 1e12 / GT_PER_MM_SLE
        total_mass = float(assemble(h * dx)) * RHO_I_SI / 1e12

        results.append((t_yr, vaf, total_mass))

        if k % output_interval == 0 or k == 1:
            PETSc.Sys.Print(
                f"  t={t_yr:.1f}  VAF={vaf:.4f} mm SLE  "
                f"mass={total_mass:.1f} Gt  [{t_elapsed:.1f}s]"
            )

        if k % checkpoint_interval == 0:
            chk_fn = os.path.join(
                RESULTS_DIR, f"{experiment_name}_{lc}_t{t_yr:.0f}.h5"
            )
            with fd.CheckpointFile(chk_fn, "w") as chk:
                chk.save_mesh(mesh)
                chk.save_function(h, name="thickness")
                chk.save_function(z.subfunctions[0], name="velocity")
            PETSc.Sys.Print(f"    [checkpoint: {chk_fn}]")

    PETSc.Sys.Print(f"\n{experiment_name} simulation complete.")

    chk_fn = os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_final.h5")
    with fd.CheckpointFile(chk_fn, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(h, name="thickness")
        chk.save_function(s, name="surface")
        chk.save_function(z.subfunctions[0], name="velocity")
    PETSc.Sys.Print(f"Saved: {chk_fn}")

    csv_fn = os.path.join(RESULTS_DIR, f"{experiment_name}_{lc}_timeseries.csv")
    with open(csv_fn, "w") as f:
        f.write("year,vaf_mm_sle,mass_gt\n")
        for t, vaf, mass in results:
            f.write(f"{t:.1f},{vaf:.6f},{mass:.2f}\n")
    PETSc.Sys.Print(f"Saved: {csv_fn}")

    return results
