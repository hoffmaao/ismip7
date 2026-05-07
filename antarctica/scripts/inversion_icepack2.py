#!/usr/bin/env python3
"""
Joint inversion for basal friction and fluidity using icepack2 + tlm_adjoint.

Follows the Kangerd demo pattern from Shapero's dual-problems repo:
- 3-field Z = V * Sigma * T (no DG thickness in mixed space)
- Regularized Jacobian: J = J_r + α * J_1
- snes_divergence_tolerance = -1
- Sliding coefficient includes exp(m*theta)

Usage:
    python scripts/inversion_icepack2.py
    mpiexec -n 16 python scripts/inversion_icepack2.py
"""

import numpy as np
import os, glob, json
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
    ds,
    dS,
    split,
    assemble,
    Mesh,
    FunctionSpace,
    VectorFunctionSpace,
    TensorFunctionSpace,
    FiniteElement,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    COMM_WORLD,
    exp,
)
from tlm_adjoint.firedrake import (
    reset_manager,
    start_manager,
    stop_manager,
    clear_caches,
    compute_gradient,
    Functional,
    EquationSolver,
)
from firedrake.petsc import PETSc
from scipy.optimize import minimize as scipy_minimize

import rasterio, icepack
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
    glen_flow_law as n_glen_val,
    weertman_sliding_law as m_weertman,
)
import colorcet as cc
import matplotlib.pyplot as plt

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
FIG_DIR = os.path.join(_ROOT, "figs")

lc = int(os.environ.get("ISMIP7_LC", "8000"))
lc_coarse = int(os.environ.get("ISMIP7_LC_COARSE", str(lc * 10)))

# Regularization
GAMMA_THETA = 1.0
GAMMA_PHI = 1.0
L_REG = 7.5e3


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

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

    use_calving_terminus = os.environ.get("ISMIP7_NO_CALVING_TERMINUS") is None

    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Sigma = TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T

    # ── Load Data ──
    PETSc.Sys.Print("Loading data...")
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    b = icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:bed"), Q)
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "10.0"))
    H = Function(Q).interpolate(
        max_value(
            icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q),
            Constant(h_clamp),
        )
    )
    PETSc.Sys.Print(f"  H clamp: {h_clamp} m")
    rho_ratio = Constant(917.0 / 1024.0)

    rho_ratio = Constant(917.0 / 1024.0)
    s = Function(Q).interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))

    vel_fn = find_file(os.path.join(DATA_DIR, "velocity"), "*.nc")
    u_obs = icepack.interpolate(
        (rasterio.open(f"netcdf:{vel_fn}:VX"), rasterio.open(f"netcdf:{vel_fn}:VY")),
        V,
        fillvalue=0.0,
    )

    calving_ids = tuple(bnd_ids["calving"])

    # ── Rheology ──
    A0 = Constant(icepack.rate_factor(Constant(260.0)))
    n = Constant(n_glen_val)
    m = Constant(n_glen_val)
    tau_c = Constant(0.1)

    u_speed = Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    )
    u_c = Constant(float(u_speed.dat.data_ro.mean()))
    PETSc.Sys.Print(f"  tau_c={float(tau_c):.3f} MPa, u_c={float(u_c):.1f} m/yr")

    sparams = {
        "snes_type": "newtonls",
        "snes_max_it": 200,
        "snes_linesearch_type": "nleqerr",
        "snes_divergence_tolerance": -1,
        "snes_stol": 0.0,
        "ksp_type": "gmres",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 400,  # working memory increase
        "mat_mumps_icntl_24": 1,  # detect null pivots
        "mat_mumps_cntl_3": 1e-12,  # null pivot threshold
    }
    fc_params = {"quadrature_degree": 4}

    # ── Build form (Kangerd pattern: controls baked into sliding coefficient) ──
    z = Function(Z)
    z.sub(0).interpolate(Constant(0.1) * u_obs)

    theta = Function(Q, name="theta")  # log friction adjustment
    phi = Function(Q, name="phi")  # log fluidity adjustment

    # Warm-start theta/phi from a previous checkpoint if available
    warm_chk = os.environ.get("ISMIP7_WARM_START")
    if warm_chk:
        PETSc.Sys.Print(f"  Loading warm start from {warm_chk}")
        with fd.CheckpointFile(warm_chk, "r") as chk:
            chk_mesh = chk.load_mesh()
            theta_ws = chk.load_function(chk_mesh, name="log_friction")
            phi_ws = chk.load_function(chk_mesh, name="log_fluidity")
        theta.dat.data[:] = theta_ws.dat.data_ro
        phi.dat.data[:] = phi_ws.dat.data_ro
        PETSc.Sys.Print(
            f"    theta: [{theta.dat.data_ro.min():.3f}, {theta.dat.data_ro.max():.3f}]"
        )
        PETSc.Sys.Print(
            f"    phi:   [{phi.dat.data_ro.min():.3f}, {phi.dat.data_ro.max():.3f}]"
        )

    u_s, M_s, tau_s = split(z)

    fields = {
        "velocity": u_s,
        "membrane_stress": M_s,
        "basal_stress": tau_s,
        "thickness": H,
        "surface": s,
    }

    # Sliding coefficient: Weertman with smooth Heaviside grounding mask
    # K = He * u_c / tau_c^n  (He = 0 floating, 1 grounded)
    # In the dual form friction_power: Ψ*(τ) = |τ|^{1+n} / ((1+n)*K)
    # He → 0 makes K → 0 which makes Ψ* → ∞, penalizing any basal
    # stress on floating ice. Add a small floor to keep K > 0.
    # Combined Ua approach: He multiplies τ in momentum_balance AND
    # K scales with 1/He so friction_power makes τ cheap on floating ice.
    # K large → τ unconstrained (floating), K small → τ penalized (grounded)
    # Together with floating=He in momentum_balance, this gives zero
    # effective friction on floating ice.
    # Smooth Heaviside sliding coefficient.
    # In the dual form: friction_power = K/(m+1)*|τ|^{m+1}
    #   Large K → τ penalized → zero friction (floating)
    #   Small K → τ allowed → normal friction (grounded)
    # He ≈ 1 grounded, He ≈ 0 floating. We use 1/(He + floor) to make
    # K large on floating ice. Floor = 0.0001 → K is 10000x larger on
    # shelves than grounded → effectively zero shelf friction.
    phi_eff = Function(Q).interpolate(
        max_value(
            Constant(1.0)
            - rho_W
            * g
            * max_value(Constant(0.0), -b)
            / (rho_I * g * max_value(H, Constant(1.0))),
            Constant(0.01),
        )
    )
    K_base = u_c / (phi_eff * tau_c) ** n

    rheology = {
        "flow_law_exponent": n,
        "flow_law_coefficient": A0 * exp(phi),
        "sliding_exponent": n,
        "sliding_coefficient": K_base * exp(-n * theta),
    }

    L = (
        model.minimization.viscous_power(**fields, **rheology)
        + model.minimization.friction_power(**fields, **rheology)
        + model.minimization.momentum_balance(**fields)
    )
    if use_calving_terminus:
        L += model.minimization.calving_terminus(**fields, outflow_ids=calving_ids)
        PETSc.Sys.Print("  Using calving_terminus BC")
    else:
        PETSc.Sys.Print("  NO calving_terminus BC (buffered mesh, h=0 at front)")
    F = derivative(L, z)

    # ── Warm start ──
    stop_manager()
    prob = NonlinearVariationalProblem(F, z, form_compiler_parameters=fc_params)
    slvr = NonlinearVariationalSolver(prob, solver_parameters=sparams)
    # Always use continuation — single solve at n=3 can fail with
    # checkpoint parameters that create ill-conditioned systems
    PETSc.Sys.Print("Warm start (continuation n=1→3)...")
    for exponent in np.linspace(1.0, n_glen_val, 5):
        n.assign(exponent)
        slvr.solve()
    PETSc.Sys.Print("  Done")

    u_init = z.subfunctions[0]
    u_mag = Function(Q).interpolate(sqrt(inner(u_init, u_init)))
    PETSc.Sys.Print(f"  u_max = {float(u_mag.dat.data_ro.max()):.0f} m/yr")

    # ── Forward function for tlm_adjoint ──
    area_val = assemble(Constant(1.0) * dx(mesh))

    def forward(theta_ctrl, phi_ctrl):
        clear_caches()
        rheology_ctrl = {
            "flow_law_exponent": n,
            "flow_law_coefficient": A0 * exp(phi_ctrl),
            "sliding_exponent": n,
            "sliding_coefficient": K_base * exp(-n * theta_ctrl),
        }
        L_ctrl = (
            model.minimization.viscous_power(**fields, **rheology_ctrl)
            + model.minimization.friction_power(**fields, **rheology_ctrl)
            + model.minimization.momentum_balance(**fields)
        )
        if use_calving_terminus:
            L_ctrl += model.minimization.calving_terminus(
                **fields, outflow_ids=calving_ids
            )
        F_ctrl = derivative(L_ctrl, z)
        # Continuation inside annotation for robustness
        for exponent in np.linspace(1.0, n_glen_val, 5):
            n.assign(exponent)
            EquationSolver(
                F_ctrl == 0,
                z,
                solver_parameters=sparams,
                form_compiler_parameters=fc_params,
            ).solve()

        u_sol, _, _ = split(z)
        J = Functional(name="J")
        J.assign(
            0.5
            / area_val
            * ((u_sol[0] - u_obs[0]) ** 2 + (u_sol[1] - u_obs[1]) ** 2)
            * dx
        )
        return J

    # ── MPI helpers ──
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

    # ── L-BFGS-B Inversion ──
    max_iter = int(os.environ.get("ISMIP7_MAXITER", "500"))
    PETSc.Sys.Print(f"\nStarting L-BFGS-B inversion (theta + phi)...")
    PETSc.Sys.Print(f"  maxiter={max_iter}, nranks={COMM_WORLD.size}")

    global_ndof = len(func_to_global(theta))
    z_backup = z.copy(deepcopy=True)
    last_good_obj = [np.inf]
    iteration_count = [0]

    gamma_theta = Constant(GAMMA_THETA)
    gamma_phi = Constant(GAMMA_PHI)
    L_reg_c = Constant(L_REG)

    def objective_and_gradient(x_vec):
        t_iter = perf_counter()
        global_to_func(x_vec[:global_ndof], theta)
        global_to_func(x_vec[global_ndof:], phi)

        t_fwd = perf_counter()
        reset_manager()
        start_manager()
        try:
            J = forward(theta, phi)
        except fd.ConvergenceError:
            stop_manager()
            z.assign(z_backup)
            PETSc.Sys.Print("  [!] Forward solve failed, returning large objective")
            return last_good_obj[0] * 10, np.zeros(2 * global_ndof)
        stop_manager()
        J_val = float(J)
        t_fwd = perf_counter() - t_fwd

        z_backup.assign(z)
        last_good_obj[0] = J_val

        t_adj = perf_counter()
        dJ_dtheta, dJ_dphi = compute_gradient(J, [theta, phi])
        t_adj = perf_counter() - t_adj

        reg_theta = float(
            assemble(
                0.5
                / area_val
                * gamma_theta
                * L_reg_c**2
                * inner(grad(theta), grad(theta))
                * dx
            )
        )
        reg_phi = float(
            assemble(
                0.5
                / area_val
                * gamma_phi
                * L_reg_c**2
                * inner(grad(phi), grad(phi))
                * dx
            )
        )
        dR_theta = assemble(
            1.0
            / area_val
            * gamma_theta
            * L_reg_c**2
            * inner(grad(theta), grad(fd.TestFunction(Q)))
            * dx
        )
        dR_phi = assemble(
            1.0
            / area_val
            * gamma_phi
            * L_reg_c**2
            * inner(grad(phi), grad(fd.TestFunction(Q)))
            * dx
        )

        g_theta = func_to_global(dJ_dtheta) + func_to_global(dR_theta)
        g_phi = func_to_global(dJ_dphi) + func_to_global(dR_phi)

        total = J_val + reg_theta + reg_phi
        total_grad = np.concatenate([g_theta, g_phi])

        t_iter = perf_counter() - t_iter
        iteration_count[0] += 1
        PETSc.Sys.Print(
            f"  iter {iteration_count[0]:3d}: "
            f"misfit={J_val:.6e} reg_θ={reg_theta:.4e} reg_φ={reg_phi:.4e} "
            f"total={total:.6e} |grad|={np.linalg.norm(total_grad):.4e} "
            f"[fwd={t_fwd:.1f}s adj={t_adj:.1f}s total={t_iter:.1f}s]"
        )

        # Periodic checkpoint every 20 iterations
        if iteration_count[0] % 20 == 0:
            chk_fn = os.path.join(MESH_DIR, f"inversion_icepack2_{lc}.h5")
            with fd.CheckpointFile(chk_fn, "w") as chk:
                chk.save_mesh(mesh)
                chk.save_function(theta, name="log_friction")
                chk.save_function(phi, name="log_fluidity")
                chk.save_function(u_obs, name="velocity_obs")
                chk.save_function(H, name="thickness")
                chk.save_function(b, name="bed")
                chk.save_function(s, name="surface")
            PETSc.Sys.Print(f"    [checkpoint saved: iter {iteration_count[0]}]")

        return total, total_grad

    x0 = np.concatenate([func_to_global(theta), func_to_global(phi)])
    result = scipy_minimize(
        objective_and_gradient,
        x0,
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": max_iter, "ftol": 0, "gtol": 0, "disp": False},
    )

    PETSc.Sys.Print(f"\nOptimization finished: {result.message}")
    PETSc.Sys.Print(f"  {result.nit} iterations, {result.nfev} function evaluations")
    global_to_func(result.x[:global_ndof], theta)
    global_to_func(result.x[global_ndof:], phi)
    PETSc.Sys.Print(
        f"  theta range: [{float(theta.dat.data_ro.min()):.3f}, {float(theta.dat.data_ro.max()):.3f}]"
    )
    PETSc.Sys.Print(
        f"  phi range:   [{float(phi.dat.data_ro.min()):.3f}, {float(phi.dat.data_ro.max()):.3f}]"
    )

    # ── Save MAP immediately ──
    chk_fn = os.path.join(MESH_DIR, f"inversion_icepack2_{lc}.h5")
    with fd.CheckpointFile(chk_fn, "w") as chk:
        chk.save_mesh(mesh)
        chk.save_function(theta, name="log_friction")
        chk.save_function(phi, name="log_fluidity")
        chk.save_function(u_obs, name="velocity_obs")
        chk.save_function(H, name="thickness")
        chk.save_function(b, name="bed")
        chk.save_function(s, name="surface")
    PETSc.Sys.Print(f"Saved MAP: {chk_fn}")

    # ── Final forward solve ──
    PETSc.Sys.Print("\nFinal forward solve...")
    stop_manager()
    try:
        slvr.solve()
    except fd.ConvergenceError:
        PETSc.Sys.Print("  Final solve failed, using last optimization state")

    u_sol = z.subfunctions[0]
    u_sol_mag = Function(Q).interpolate(sqrt(u_sol[0] ** 2 + u_sol[1] ** 2))
    misfit = float(
        assemble(
            0.5
            / area_val
            * ((u_sol[0] - u_obs[0]) ** 2 + (u_sol[1] - u_obs[1]) ** 2)
            * dx
        )
    )
    PETSc.Sys.Print(f"  Final misfit: {misfit:.6e}")

    # Update checkpoint with velocity
    with fd.CheckpointFile(chk_fn, "a") as chk:
        chk.save_function(u_sol, name="velocity")
    PETSc.Sys.Print(f"Saved velocity: {chk_fn}")

    # ── Plot ──
    PETSc.Sys.Print("Plotting...")
    u_speed = Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    )
    coords = mesh.coordinates.dat.data_ro
    xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
    ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
    pad = 0.02 * max(xmax - xmin, ymax - ymin)

    pw, ph, bot = 0.25, 0.75, 0.10
    gap1, gap2 = 0.03, 0.10
    cb_w, cb_gap = 0.012, 0.015
    cb_h = ph * 0.5
    cb_bot = bot + (ph - cb_h) / 2
    x0 = 0.04
    x1 = x0 + pw + gap1
    x_cb1 = x1 + pw + cb_gap
    x2 = x_cb1 + cb_w + gap2
    x_cb2 = x2 + pw + cb_gap

    fig = plt.figure(figsize=(22, 7))
    ax0 = fig.add_axes([x0, bot, pw, ph])
    ax1 = fig.add_axes([x1, bot, pw, ph])
    cax1 = fig.add_axes([x_cb1, cb_bot, cb_w, cb_h])
    ax2 = fig.add_axes([x2, bot, pw, ph])
    cax2 = fig.add_axes([x_cb2, cb_bot, cb_w, cb_h])
    for ax in [ax0, ax1, ax2]:
        ax.set_aspect("equal")
        ax.set_xlim(xmin - pad, xmax + pad)
        ax.set_ylim(ymin - pad, ymax + pad)
        fd.triplot(
            mesh,
            axes=ax,
            interior_kw={"linewidth": 0.1, "alpha": 0.3, "color": "k"},
            boundary_kw={"linewidth": 1.0, "color": "k"},
        )
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])

    fd.tripcolor(u_speed, vmin=0, vmax=1000, axes=ax0, cmap=cc.cm.CET_L19)
    ax0.set_title("Observed speed")
    cs = fd.tripcolor(u_sol_mag, vmin=0, vmax=1000, axes=ax1, cmap=cc.cm.CET_L19)
    ax1.set_title("Inverted speed (icepack2)")
    diff = Function(Q).interpolate(u_speed - u_sol_mag)
    cd = fd.tripcolor(diff, vmin=-500, vmax=500, axes=ax2, cmap=cc.cm.CET_CBTD1)
    ax2.set_title("Observed - Modeled")
    fig.colorbar(cs, cax=cax1, label="m/yr")
    fig.colorbar(cd, cax=cax2, label="m/yr")

    out_fn = os.path.join(FIG_DIR, f"inversion_icepack2_{lc}.png")
    fig.savefig(out_fn, dpi=200, bbox_inches="tight")
    PETSc.Sys.Print(f"Saved: {out_fn}")


if __name__ == "__main__":
    main()
