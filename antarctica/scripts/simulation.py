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

# SI densities for diagnostics (icepack2 constants are in MPa-m-yr units)
_RHO_I_SI = 917.0
_RHO_W_SI = 1024.0

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")
RESULTS_DIR = os.path.join(_ROOT, "results")

# Repo root on the path for the shared dual-friction operator.
sys.path.insert(0, os.path.dirname(_ROOT))

lc = int(os.environ.get("ISMIP7_LC", "2500"))
lc_coarse = int(os.environ.get("ISMIP7_LC_COARSE", "64000"))


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

    # Friction law: "budd" (power-law dual, default) or "regularized_coulomb"
    # (RC residual driven by the inversion_icepack2_rc_<lc>.h5 MAP:
    # grounded-only theta, exact-zero shelf drag via the Coulomb cap c0*N).
    # Must mirror inversion_icepack2.py so theta/phi keep their meaning.
    friction = os.environ.get("ISMIP7_FRICTION", "budd")
    use_rc = friction == "regularized_coulomb"
    map_tag = "_rc" if use_rc else ""

    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Sigma = TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T

    PETSc.Sys.Print("Loading data...")
    bm_fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    b = icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:bed"), Q)
    # Two clamps:
    #   - h_clamp_init: floor on the *initial* thickness from BedMachine so
    #     the initial diagnostic solve sees a well-posed problem everywhere
    #     (default 10 m, matches the historical behavior).
    #   - h_clamp: floor enforced in the advection step. Default 0 so melt
    #     can drive cells to zero at the calving front; the composite
    #     rheology keeps the diagnostic SNES nonsingular at h=0.
    # In RC mode the MAP was inverted against the TRUE h=0 BedMachine
    # geometry (h_visc_floor handles the ice-free buffer), so the initial
    # clamp defaults to 0 and the initial state never touches a floor.
    h_clamp_init = float(
        os.environ.get("ISMIP7_H_CLAMP_INIT", "0.0" if use_rc else "10.0")
    )
    h_clamp = float(os.environ.get("ISMIP7_H_CLAMP", "0.0"))
    H = Function(Q, name="thickness").interpolate(
        max_value(
            icepack.interpolate(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q),
            Constant(h_clamp_init),
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

    inv_fn = os.path.join(MESH_DIR, f"inversion_icepack2{map_tag}_{lc}.h5")
    PETSc.Sys.Print(f"Loading MAP: {inv_fn}")
    with fd.CheckpointFile(inv_fn, "r") as chk:
        chk_mesh = chk.load_mesh()
        # The dof copies below require identical node ordering between the
        # .msh mesh and the checkpoint mesh; verify via coordinates.
        _co, _cc = mesh.coordinates.dat.data_ro, chk_mesh.coordinates.dat.data_ro
        if _cc.shape != _co.shape or not np.allclose(_cc, _co):
            raise RuntimeError(
                f"{inv_fn}: checkpoint mesh nodes do not match {mesh_fn}; "
                "set ISMIP7_MESH to the mesh the inversion ran on"
            )
        theta = chk.load_function(chk_mesh, name="log_friction")
        phi = chk.load_function(chk_mesh, name="log_fluidity")
        if use_rc:
            H_chk = chk.load_function(chk_mesh, name="thickness")
            b_chk = chk.load_function(chk_mesh, name="bed")
            s_chk = chk.load_function(chk_mesh, name="surface")
            u_obs_chk = chk.load_function(chk_mesh, name="velocity_obs")
    theta_f = Function(Q, name="theta")
    theta_f.dat.data[:] = theta.dat.data_ro
    phi_f = Function(Q, name="phi")
    phi_f.dat.data[:] = phi.dat.data_ro

    # Clip the log-adjustments to a sane band. theta/phi are O(1) in a
    # converged MAP (the rc_32000 MAP maxes at 1.64), so anything beyond
    # ISMIP7_MAP_CLIP (default 6) is optimization noise from an
    # unconverged checkpoint — e.g. the rc_500 MAP (a 20-iter snapshot)
    # carries ~700 nodes with |theta|,|phi| up to 23, i.e. exp() coeffs
    # ~1e10 that are local singularities the diagnostic SNES cannot solve
    # through. Clipping removes the pathology (0.1% of nodes) without
    # touching legitimate structure. Set 0 to disable.
    map_clip = float(os.environ.get("ISMIP7_MAP_CLIP", "6.0"))
    if map_clip > 0.0:
        n_clip = 0
        for fld in (theta_f, phi_f):
            d = fld.dat.data
            n_clip += int((np.abs(d) > map_clip).sum())
            np.clip(d, -map_clip, map_clip, out=d)
        if n_clip:
            PETSc.Sys.Print(
                f"  MAP clip: bounded {n_clip} theta/phi node(s) to "
                f"|.|<={map_clip:.0f} (unconverged-checkpoint outliers)"
            )

    if use_rc:
        # Use the exact geometry/velocity the RC inversion saw, so the
        # Weertman anchor C_w0 (and hence the meaning of theta) is
        # reproduced bit-for-bit. h_clamp_init (default 0) may re-floor.
        H.dat.data[:] = H_chk.dat.data_ro
        b.dat.data[:] = b_chk.dat.data_ro
        s.dat.data[:] = s_chk.dat.data_ro
        u_obs.dat.data[:] = u_obs_chk.dat.data_ro
        if h_clamp_init > 0.0:
            H.interpolate(max_value(H, Constant(h_clamp_init)))
            s.interpolate(max_value(b + H, (Constant(1.0) - rho_ratio) * H))

    A0 = Constant(icepack.rate_factor(Constant(260.0)))
    # Composite Goldsby-Kohlstedt-style exponents (must match the inversion
    # that produced the MAP file we load above).
    n_flow_val = float(os.environ.get("ISMIP7_N_FLOW", "4.0"))
    m_slide_val = float(os.environ.get("ISMIP7_M_SLIDE", "3.0"))
    a4_factor = float(os.environ.get("ISMIP7_A4_FACTOR", "10.0"))
    n_flow = Constant(n_flow_val)
    m_slide = Constant(m_slide_val)
    tau_c = Constant(0.1)
    u_c = Constant(float(Function(Q).interpolate(
        max_value(sqrt(u_obs[0] ** 2 + u_obs[1] ** 2), Constant(1.0))
    ).dat.data_ro.mean()))

    # Phi_eff (effective-pressure fraction). Uses a small floor on H so it
    # is well-defined where the original BedMachine thickness is 0.
    _H_FLOOR_PHI = Constant(1.0)
    phi_eff = Function(Q).interpolate(
        max_value(
            Constant(1.0)
            - rho_W * g * max_value(Constant(0.0), -b)
              / (rho_I * g * max_value(H, _H_FLOOR_PHI)),
            Constant(0.01),
        )
    )
    A4_base = A0 * Constant(a4_factor)
    A_map = A4_base * exp(phi_f)
    K_base = u_c / (phi_eff * tau_c) ** m_slide
    K_map = K_base * exp(-m_slide * theta_f)

    # Composite rheology: Goldsby-Kohlstedt dislocation creep (n=n_flow=4)
    # + α · linear regularizer (n=1) with a CONSTANT reference thickness
    # H_ref. This pins M where h → 0 so the SNES Jacobian stays
    # nonsingular at the calving front and h is allowed to reach zero.
    # Mirrors icepack2 dome_test.py.
    # RC MAP was inverted under alpha=1e-2 (SNES robustness across many
    # forward solves); keep the forward diagnostic consistent with it.
    alpha_reg = Constant(float(
        os.environ.get("ISMIP7_COMPOSITE_ALPHA", "1e-2" if use_rc else "1e-4")
    ))
    H_ref = Constant(float(os.environ.get("ISMIP7_H_REF", "100.0")))
    A_linear = A_map * tau_c ** (n_flow_val - 1)   # linearized at tau_c
    K_linear = u_c / (phi_eff * tau_c) * exp(-theta_f)  # linearized at tau_c

    C_w0 = None
    if use_rc:
        from icepack2_tools.dual_friction import (
            build_rc_residual, weertman_anchor,
        )
        c0_rc = float(os.environ.get("ISMIP7_RC_C0", "0.5"))
        rc_hvisc_floor = float(os.environ.get("ISMIP7_RC_HVISC_FLOOR", "10.0"))
        rc_cw0_floor = float(os.environ.get("ISMIP7_RC_CW0_FLOOR", "0.0"))
        # Minimum Coulomb yield [MPa] -> a small, nonzero drag even where
        # N=0 (floating). RC's exact-zero shelf drag leaves floating ice
        # with no friction to damp numerical noise in the driving stress,
        # so the diagnostic↔transport coupling is unstable (velocity blows
        # up ~2x/step); a few kPa restores the damping Budd gets for free
        # from its phi_eff floor. 0 = bit-exact-zero shelves (unstable in
        # prognostic use). See the multi-step blow-up diagnosis, Jul 2026.
        rc_eps_tauc = float(os.environ.get("ISMIP7_RC_EPS_TAUC", "0.0"))
        # Anchor from the inversion-time geometry/velocity, BEFORE any
        # restart overwrites s: theta is a log-adjustment on this C_w0.
        C_w0 = weertman_anchor(H, s, u_obs, m_slide_val, Q)
        PETSc.Sys.Print(
            f"  Friction: regularized Coulomb (c0={c0_rc}, "
            f"h_visc_floor={rc_hvisc_floor:.0f}m, cw0_floor={rc_cw0_floor:.1e}, "
            f"eps_tauc={rc_eps_tauc:.1e} MPa, alpha={float(alpha_reg):.1e})"
        )

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
    # Optional SNES/KSP convergence monitoring (ISMIP7_SNES_MONITOR=1).
    # ISMIP7_SNES_LOG routes the output to a file (per run, so concurrent
    # debug runs don't interleave); otherwise it goes to stdout.
    if os.environ.get("ISMIP7_SNES_MONITOR"):
        _snes_log = os.environ.get("ISMIP7_SNES_LOG")
        _viewer = f"ascii:{_snes_log}" if _snes_log else None
        sparams.update({
            "snes_monitor": _viewer,
            "snes_converged_reason": _viewer,
            "snes_linesearch_monitor": _viewer,
            "ksp_converged_reason": _viewer,
        })
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
        h.interpolate(max_value(h, Constant(h_clamp)))  # don't restart below the floor
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
    rheo_glen = {
        "flow_law_exponent": n_flow,
        "flow_law_coefficient": A_map,
        "sliding_exponent": m_slide,
        "sliding_coefficient": K_map,
    }
    rheo_linear = {
        "flow_law_exponent": Constant(1.0),
        "flow_law_coefficient": A_linear,
        "sliding_exponent": Constant(1.0),
        "sliding_coefficient": K_linear,
    }
    # Composite-rheology fields: use the constant reference thickness for the
    # linear regularization terms so they stay positive-definite at h=0.
    fields_reg = dict(fields)
    fields_reg["thickness"] = H_ref

    if use_rc:
        # Residual closure on the LIVE prognostic fields (h, s): the
        # grounded gate, Coulomb cap, and driving stress all track the
        # evolving geometry, so the grounding line migrates freely.
        F = build_rc_residual(
            z, theta_f, phi_f, H=h, s=s, b=b, C_w0=C_w0,
            A4_base=A4_base, n_flow=n_flow, n_flow_val=n_flow_val,
            m_slide=m_slide_val, tau_c=tau_c, alpha=alpha_reg, H_ref=H_ref,
            c0=c0_rc, eps_tauc=rc_eps_tauc,
            c_w0_floor=rc_cw0_floor, h_visc_floor=rc_hvisc_floor,
            calving_ids=calving_ids if use_calving_terminus else None,
        )
    else:
        L = (
            model.minimization.viscous_power(**fields, **rheo_glen)
            + alpha_reg * model.minimization.viscous_power(**fields_reg, **rheo_linear)
            + model.minimization.friction_power(**fields, **rheo_glen)
            + alpha_reg * model.minimization.friction_power(**fields, **rheo_linear)
            + model.minimization.momentum_balance(**fields)
        )
        if use_calving_terminus:
            L += model.minimization.calving_terminus(**fields, outflow_ids=calving_ids)
        F = derivative(L, z)

    prob = NonlinearVariationalProblem(
        F, z, form_compiler_parameters=fc_params
    )
    slvr = NonlinearVariationalSolver(prob, solver_parameters=sparams)

    # Adaptive n/m continuation for the cold-start diagnostic solve. On a
    # fine mesh with a rough (mid-optimization) MAP the n=1→n_flow_val jump
    # can outrun Newton (DIVERGED_MAX_IT); restore the initial guess and
    # re-ramp with more, smaller steps rather than crashing. Escalates
    # ISMIP7_CONTINUATION_STEPS (default 8) → 2× → 4×.
    base_steps = int(os.environ.get("ISMIP7_CONTINUATION_STEPS", "8"))
    z_init = z.copy(deepcopy=True)
    PETSc.Sys.Print(
        f"Initial diagnostic solve (continuation n_flow 1→{n_flow_val:.1f}, "
        f"m_slide 1→{m_slide_val:.1f})..."
    )
    for attempt, steps in enumerate((base_steps, 2 * base_steps, 4 * base_steps)):
        try:
            for t in np.linspace(0.0, 1.0, steps):
                n_flow.assign(1.0 + t * (n_flow_val - 1.0))
                m_slide.assign(1.0 + t * (m_slide_val - 1.0))
                slvr.solve()
            PETSc.Sys.Print(f"  Done ({steps} continuation steps)")
            break
        except fd.ConvergenceError:
            if attempt == 2:
                PETSc.Sys.Print(
                    f"  Continuation diverged at {steps} steps — giving up."
                )
                raise
            PETSc.Sys.Print(
                f"  Continuation diverged at {steps} steps; "
                f"restarting with {2 * steps}..."
            )
            z.assign(z_init)
            n_flow.assign(1.0)
            m_slide.assign(1.0)

    u0 = z.subfunctions[0]
    area = assemble(Constant(1.0) * dx(mesh))
    misfit0 = float(assemble(
        0.5 / area * ((u0[0] - u_obs[0]) ** 2 + (u0[1] - u_obs[1]) ** 2) * dx
    ))
    PETSc.Sys.Print(f"  Initial velocity misfit vs obs: {misfit0:.6e}")

    accum = Function(Q, name="accumulation").assign(0.0)
    ocean_melt = Function(Q, name="ocean_melt").assign(0.0)

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
        "n_flow": n_flow,
        "n_flow_val": n_flow_val,
        "m_slide": m_slide,
        "m_slide_val": m_slide_val,
        "accum": accum,
        "ocean_melt": ocean_melt,
        "phi_eff": phi_eff,
        "rho_ratio": rho_ratio,
        "h_clamp": h_clamp,
        "calving_ids": calving_ids,
        "u_obs": u_obs,
        "friction": friction,
        # t=0 (BedMachine/inversion) thickness — NOT overwritten by
        # restart_from, so the fixed calving front stays anchored to the
        # observed extent across restarts.
        "H_init": H,
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
    n_flow = ctx["n_flow"]
    n_flow_val = ctx["n_flow_val"]
    m_slide = ctx["m_slide"]
    m_slide_val = ctx["m_slide_val"]
    accum = ctx["accum"]
    ocean_melt = ctx["ocean_melt"]
    phi_eff = ctx["phi_eff"]
    rho_ratio = ctx["rho_ratio"]
    h_clamp = ctx["h_clamp"]

    # round, don't truncate: int((2300-2015)/0.1) = 2849 loses the last step
    nsteps = int(round((t_end - t_start) / dt))
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

    rho_gt = _RHO_I_SI / 1e12  # m^3 ice -> Gt

    # Fixed calving front (ISMIP7_FIXED_FRONT=1): cells that are ice-free
    # in the initial state may not accumulate ice; whatever flows into
    # them is removed each step and tallied as calving flux. Without this
    # a buffered mesh has NO calving sink (~1300 Gt/yr in reality) and the
    # sheet must gain mass. Only meaningful when the initial state is the
    # true BedMachine geometry (RC mode / h_clamp_init=0) — with a clamped
    # initial state every cell has ice and the mask is empty.
    fixed_front = os.environ.get("ISMIP7_FIXED_FRONT") is not None
    front_hmin = float(os.environ.get("ISMIP7_FRONT_HMIN", "1.0"))
    beyond_front = None
    cell_area = assemble(fd.TestFunction(Q_dg) * dx).dat.data_ro.copy()
    if fixed_front:
        # Mask from the t=0 observed extent (ctx["H_init"]), not the
        # current h: a restarted run must not re-mask cells that
        # legitimately retreated mid-run inside the observed extent.
        h_dg.project(ctx.get("H_init", h))
        beyond_front = h_dg.dat.data_ro < front_hmin
        n_beyond = mesh.comm.allreduce(int(beyond_front.sum()))
        PETSc.Sys.Print(
            f"  Fixed calving front: {n_beyond} initially ice-free cells "
            f"masked (h < {front_hmin} m)"
        )

    mass_prev = float(assemble(h * dx)) * rho_gt

    # ISMIP7_LEGACY_TRANSPORT=1 restores the pre-Jul-2026 scheme: the
    # -h*div(u*phi) volume term (non-conservative for DG0: it adds
    # spurious h*div(u) mass at thickness jumps) and the L2 DG0->CG1
    # projection (overshoots negative at fronts; the h floor then
    # injects mass). The default is the exactly-conservative FV form
    # plus a lumped-mass projection (convex combination of adjacent
    # cell values: bounded and integral-preserving).
    legacy_transport = os.environ.get("ISMIP7_LEGACY_TRANSPORT") is not None
    if legacy_transport:
        PETSc.Sys.Print("  LEGACY transport: non-conservative volume term + L2 projection")
    m_lump = assemble(fd.TestFunction(Q) * dx)
    proj_rhs = fd.Cofunction(Q.dual())

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
            try:
                for t in np.linspace(0.0, 1.0, 10):
                    n_flow.assign(1.0 + t * (n_flow_val - 1.0))
                    m_slide.assign(1.0 + t * (m_slide_val - 1.0))
                    slvr.solve()
            except fd.ConvergenceError:
                PETSc.Sys.Print(f"  Step {k}: continuation failed, saving and stopping")
                break

        u_vel = z.subfunctions[0]

        # DG0 upwind implicit Euler for thickness
        h_dg.project(h)
        h_dg_old.assign(h_dg)

        un = fd.dot(u_vel, n_facet)
        un_plus = (un + abs(un)) / 2

        F_prog = (
            (h_dg_trial - h_dg_old) / dt_c * phi_dg * dx
            + (un_plus("+") * h_dg_trial("+") - un_plus("-") * h_dg_trial("-"))
            * fd.jump(phi_dg)
            * dS
            + un_plus * h_dg_trial * phi_dg * ds
            - (accum - ocean_melt) * phi_dg * dx
        )
        if legacy_transport:
            F_prog += -h_dg_trial * fd.div(u_vel * phi_dg) * dx
        fd.solve(
            fd.lhs(F_prog) == fd.rhs(F_prog),
            h_dg,
            solver_parameters={
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",
            },
        )

        # Mass budget: attribute this step's dM to its sources so a drift
        # is diagnosable (SMB - melt - boundary outflux - calving) and any
        # non-conservative floor/projection mass shows up explicitly.
        smb_rate = float(assemble(accum * dx)) * rho_gt          # Gt/yr
        melt_rate = float(assemble(ocean_melt * dx)) * rho_gt    # Gt/yr
        out_rate = float(assemble(un_plus * h_dg * ds)) * rho_gt # Gt/yr
        m1 = float(assemble(h_dg * dx)) * rho_gt

        h_dg.interpolate(max_value(h_dg, Constant(h_clamp)))
        m2 = float(assemble(h_dg * dx)) * rho_gt
        clamp_gt = m2 - m1                                       # Gt added by DG floor

        calv_gt = 0.0
        if fixed_front:
            data = h_dg.dat.data
            calv_gt = mesh.comm.allreduce(
                float((data[beyond_front] * cell_area[beyond_front]).sum())
            ) * rho_gt
            data[beyond_front] = 0.0

        if legacy_transport:
            h.project(h_dg)
        else:
            # Lumped-mass projection: h_i = ∫phi_i h_dg / ∫phi_i.
            assemble(fd.TestFunction(Q) * h_dg * dx, tensor=proj_rhs)
            h.dat.data[:] = proj_rhs.dat.data_ro / m_lump.dat.data_ro
        m3 = m2 - calv_gt                                        # ∫h preserved by projection
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
        vaf = float(assemble(haf * dx)) * _RHO_I_SI / 1e12 / 362.5
        total_mass = float(assemble(h * dx)) * _RHO_I_SI / 1e12

        clamp_cg_gt = total_mass - m3        # Gt added by the CG floor after projection
        dm = total_mass - mass_prev
        resid_gt = dm - (
            (smb_rate - melt_rate - out_rate) * dt
            + clamp_gt + clamp_cg_gt - calv_gt
        )
        mass_prev = total_mass

        results.append((t_yr, vaf, total_mass, smb_rate, melt_rate,
                        out_rate, calv_gt, clamp_gt + clamp_cg_gt, resid_gt))

        if k % output_interval == 0 or k == 1:
            PETSc.Sys.Print(
                f"  t={t_yr:.1f}  VAF={vaf:.4f} mm SLE  "
                f"mass={total_mass:.1f} Gt  [{t_elapsed:.1f}s]\n"
                f"      budget [Gt/yr]: SMB={smb_rate:+.0f} melt={-melt_rate:+.0f} "
                f"outflux={-out_rate:+.0f} calv={-calv_gt/dt:+.0f} "
                f"clamp={(clamp_gt+clamp_cg_gt)/dt:+.1f} "
                f"dM/dt={dm/dt:+.0f} resid={resid_gt/dt:+.2f}"
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
    if mesh.comm.rank == 0:
        with open(csv_fn, "w") as f:
            f.write("year,vaf_mm_sle,mass_gt,smb_gtyr,melt_gtyr,"
                    "outflux_gtyr,calv_gt,clamp_gt,resid_gt\n")
            for row in results:
                f.write(
                    f"{row[0]:.1f},{row[1]:.6f},{row[2]:.2f},"
                    + ",".join(f"{v:.4f}" for v in row[3:]) + "\n"
                )
    PETSc.Sys.Print(f"Saved: {csv_fn}")

    return results
