#!/usr/bin/env python3
r"""Leading eigenmodes of the prior-preconditioned Gauss-Newton Hessian of a MAP.

The Laplace approximation of Recinos et al. (2023): at the MAP
the posterior covariance is ``(H_GN + Gamma^-1)^-1`` with ``H_GN`` the
Gauss-Newton Hessian of the data term and ``Gamma`` the prior covariance.
Its departure from the prior lives in the generalised eigenpairs

    H_GN v_i = lambda_i Gamma^-1 v_i,

ordered by lambda_i: a mode with lambda_i >> 1 is constrained by the data,
one with lambda_i << 1 is known no better than the prior. Those pairs are what
this script computes and saves, and what a forward's QoI sensitivity is
projected onto (``sigma_Q^2 = sigma_prior^2 - sum_i (g.v_i)^2 lambda_i/(1+lambda_i)``).

The problem is the MAP's own: ``simulation.setup_model()`` loads the MAP
(``ISMIP7_INVERSION``) exactly as a forward does - its mesh, geometry, friction
control (sqrt, exp or log), sub-element scheme and the mixed state it was
accepted at - and the data term is the inversion's (the per-datum chi^2 on
the MEaSUREs errors with the MAP's recorded ``misfit_scale`` and log-speed
weight). The Hessian is tlm_adjoint's Gauss-Newton action ``J^T R^-1 J`` on
the recorded state, ``R^-1`` being that data term's second derivative with
respect to the velocity: positive semi-definite by construction, with no
dependence on the residual. (Differentiating the data term itself, evaluated
against the MAP's own velocity so the residual vanished, gave the same matrix
in exact arithmetic but an identically zero action in tlm_adjoint, whose
first-order adjoint is then zero and is skipped.) The prior is the MAP's
recorded bi-Laplacian (sigma, rho per control).

The tape holds one solve at the physical exponents from the converged state
(the direct forward of 1 Oct 2026); the n=1->3 ladder inside the tape is
what made the earlier Hessian indefinite. Run it under the full mixed-Jacobian
MUMPS the inversion differentiates through:

    ISMIP7_INVERSION=<map.h5> ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=full_mumps \
        mpiexec -n 4 python antarctica/scripts/run_eigendec.py [--modes 40]

Output beside the MAP's results directory: ``eigendec_<lc>[_<tag>].h5`` (the
modes as alpha/phi pairs, Gamma^-1-orthonormal), ``eigenvalues_<lc>[_<tag>].txt``
and ``eigendec_<lc>[_<tag>].json`` (the prior and likelihood the modes belong
to). The eigensolve is ARPACK on rank 0 over the global control vector;
every rank serves the collective Hessian and prior actions.

Uncertainty of the discharge (``--qoi discharge``, the default): on the same
recorded solve, the ice leaving grounded cells at the MAP, across the
grounding line into floating cells and across the mesh boundary
(``icepack2_tools.laplace.grounded_discharge_forms``), in Gt/yr, in total and
per IMBIE2 basin, and the grounding-line flux alone. Each gets its gradient
with respect to (alpha, phi) from one adjoint, its prior variance
``g . Gamma g`` before the eigensolve and its Laplace posterior variance
after it (``icepack2_tools.laplace.posterior_variance``), written to
``uncertainty_<lc>[_<tag>].json``; the prior-weighted sensitivity
``Gamma g`` of the totals goes into the eigendecomposition file. ``--check``
compares the total's adjoint gradient with a central difference of two
untaped solves along ``Gamma g``, and the Gauss-Newton action along that
direction with the norm of the velocity change two more solves give; it
stops after the checks.

``--likelihood map`` (the default) takes the Hessian of the MAP's own
objective (``map_objective_terms``); ``--likelihood velocity`` keeps the
velocity term alone, so for a MAP that also scored dH/dt (``dhdt_weight`` > 0
or ``dhdt_reweight`` > 0) its posterior is wider than the inversion's: an
upper bound. Every output records the likelihood terms and their settings.
"""
import argparse
import json
import os
import resource
import sys
import traceback
from time import perf_counter
from types import SimpleNamespace

import numpy as np
import firedrake as fd
from firedrake import (Constant, Function, TestFunction,
                       TrialFunction, assemble, dx, inner, max_value, sqrt,
                       VectorFunctionSpace)
from firedrake.petsc import PETSc
from mpi4py import MPI
from tlm_adjoint.firedrake import (CachedGaussNewton, EquationSolver, Functional,
                                   compute_gradient, reset_manager, start_manager,
                                   stop_manager)

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

import icepack  # noqa: E402
import rasterio  # noqa: E402
import simulation  # noqa: E402
from icepack2_tools.prior import (bilaplacian_coeffs,  # noqa: E402
                                  prior_operator_form)
from icepack2_tools.forcing import _basin_on_mesh, imbie2_basin_path  # noqa: E402
from icepack2_tools.laplace import (arpack_action_bound, arpack_ncv, grounded_discharge_forms,  # noqa: E402
                                    grounded_indicator, is_dg0, leading_eigenpairs, posterior_variance,
                                    tangent_linear_release_check)
from icepack2_tools.runconfig import obs_data_root  # noqa: E402
from icepack2_tools.transport import transport_form  # noqa: E402

TRANSPORT_SP = {"snes_type": "ksponly", "ksp_type": "gmres", "pc_type": "bjacobi",
                "sub_pc_type": "ilu", "ksp_rtol": 1e-10}           # the inversion's
FILTER_SP = {"snes_type": "ksponly", "ksp_type": "cg", "pc_type": "gamg",
             "ksp_rtol": 1e-12, "ksp_atol": 0.0}                   # the inversion's

RHO_GT = 917.0 / 1e12          # m3 of ice -> Gt


def rss_mb():
    """This rank's resident set in MB: the current one where /proc exists,
    the peak from getrusage elsewhere."""
    try:
        with open("/proc/self/statm") as f:
            return int(f.read().split()[1]) * resource.getpagesize() / 2 ** 20
    except OSError:
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def dual_dot(g, v):
    """The pairing of an assembled dual with a primal function, over ranks."""
    with g.dat.vec_ro as gv, v.dat.vec_ro as vv:
        return gv.dot(vv)


def discharge_quantities(ctx, u):
    """(name, UFL form in Gt/yr) for the ice leaving grounded cells at the MAP:
    the total, the grounding-line flux alone and the total per IMBIE2 basin.
    Grounded cells are those with a surface above the flotation surface, and
    every cell on a bed at or above sea level, as the ISMIP7 writer books
    ``ligroundf``, so ice flowing onto bare land is not discharge."""
    h, s, b = ctx["h"], ctx["s"], ctx["b"]
    Q0 = h.function_space()
    if not is_dg0(Q0):
        raise SystemExit("run_eigendec.py --qoi discharge reads the DG0 transport's fluxes; "
                         f"this MAP keeps its thickness in {Q0.ufl_element()}")
    grounded = grounded_indicator(Q0, s, b, ctx["rho_ratio"])
    basin = _basin_on_mesh(ctx["geom_xy"][0], ctx["geom_xy"][1], imbie2_basin_path())
    comm = Q0.mesh().comm
    present = sorted(set().union(*comm.allgather(
        set(int(v) for v in np.unique(basin[grounded.dat.data_ro > 0.5]) if v >= 0))))
    gl, edge = grounded_discharge_forms(u, h, grounded)
    out = [("discharge_total", Constant(RHO_GT) * (gl + edge)),
           ("grounding_line_flux_total", Constant(RHO_GT) * gl)]
    for k in present:
        w = Function(Q0, name=f"basin_{k:02d}")
        w.dat.data[:] = (basin == k).astype(float)
        gl_k, edge_k = grounded_discharge_forms(u, h, grounded, w)
        out.append((f"discharge_basin_{k:02d}", Constant(RHO_GT) * (gl_k + edge_k)))
    return out


def to_global(f):
    with f.dat.vec_ro as v:
        scatter, x_seq = PETSc.Scatter.toAll(v)
        scatter.scatter(v, x_seq, mode=PETSc.Scatter.Mode.FORWARD)
        out = x_seq.array.copy()
        scatter.destroy()
        x_seq.destroy()
    return out


def from_global(arr, f):
    with f.dat.vec_wo as v:
        x_seq = PETSc.Vec().createSeq(len(arr), comm=PETSc.COMM_SELF)
        x_seq.array[:] = arr
        scatter, _ = PETSc.Scatter.toAll(v)
        scatter.scatter(x_seq, v, mode=PETSc.Scatter.Mode.REVERSE)
        scatter.destroy()
        x_seq.destroy()


def map_objective_terms(ctx, meta, Q, likelihood):
    r"""The likelihood's terms beyond the velocity chi-square, rebuilt as
    ``inversion_icepack2.py`` built them from the settings the MAP records
    (``meta``) and the run's environment, for ``likelihood`` ``map``; under
    ``velocity`` there are none.

    - ``vel_weight``: the dH/dt re-weighting of the velocity chi-square,
      ``1 + a min(|dH/dt_obs| / scale, 1)`` on grounded observed cells lifted
      to the nodes (``dhdt_reweight`` a, ``dhdt_reweight_scale``), or None.
    - ``dhdt``: the dH/dt term (``dhdt_weight`` > 0), or None: one upwind
      backward-Euler step of ``ISMIP7_DHDT_DT`` years from the MAP thickness
      with the MAP velocity, its tendency masked to grounded cells Smith
      covers, the Helmholtz filter of length ``dhdt_filter_l`` normalised by
      the filtered mask, scored at ``ISMIP7_DHDT_SIGMA``, and the net term
      on the unfiltered masked tendency at ``dhdt_net_sigma`` Gt/yr.

    The step's source (SMB minus melt) and the observed dH/dt are left out:
    the Gauss-Newton Hessian does not see the residual, and the source moves
    the stepped thickness the upwind flux carries by about dt |a| / H, 1e-4
    on grounded ice. ``record`` lists the terms and their settings, written
    beside every result so modes and solves of different likelihoods are
    never combined. Collective."""
    from icepack2_tools.cellview import CellView
    from icepack2_tools.geometry import cg1_lift
    from icepack2_tools.grounding import height_above_flotation
    from icepack2_tools.obs_dhdt import grounded_cells, load_dhdt_obs

    record = {"likelihood": likelihood, "terms": ["velocity"]}
    out = SimpleNamespace(vel_weight=None, dhdt=None, record=record)
    if likelihood == "velocity":
        return out
    if likelihood != "map":
        raise ValueError(f"likelihood must be map or velocity, got {likelihood!r}")
    if "dhdt_speed_ratio" in meta:
        ratio, source = float(meta["dhdt_speed_ratio"]), "the MAP's dhdt_speed_ratio"
    else:
        ratio = float(os.environ.get("ISMIP7_DHDT_SPEED_RATIO", "1"))
        source = "ISMIP7_DHDT_SPEED_RATIO (the MAP records no dhdt_speed_ratio)"
    if ratio != 1.0:
        raise SystemExit(f"run_eigendec.py --likelihood map: {source} is {ratio:g}; the speed "
                         "weight of the velocity chi-square is not wired into the Hessian yet, "
                         "so only a ratio of 1 is accepted; use --likelihood velocity")
    w = float(meta.get("dhdt_weight", 0.0))
    a = float(meta.get("dhdt_reweight", 0.0))
    if w <= 0.0 and a <= 0.0:
        return out
    mesh = Q.mesh()
    H, b = ctx["h"], ctx["b"]
    Q_g = H.function_space()
    if not is_dg0(Q_g):
        raise SystemExit("run_eigendec.py --likelihood map rebuilds the DG0 dH/dt term; "
                         f"this MAP keeps its thickness in {Q_g.ufl_element()}")
    dobs, cov = load_dhdt_obs(Q_g)
    if a > 0.0:
        scale = float(meta.get("dhdt_reweight_scale", 1.0))
        haf = Function(Q_g).interpolate(height_above_flotation(H, b))
        k = Function(Q_g)
        k.dat.data[:] = np.where(
            (cov.dat.data_ro > 0.5) & (haf.dat.data_ro > 0.0) & (H.dat.data_ro > 1.0),
            np.minimum(np.abs(dobs.dat.data_ro) / scale, 1.0), 0.0)
        out.vel_weight = Function(Q, name="vel_weight").interpolate(
            Constant(1.0) + Constant(a) * cg1_lift(k))
        record["terms"].append("dhdt_reweight")
        record.update({"dhdt_reweight": a, "dhdt_reweight_scale": scale})
    if w > 0.0:
        grounded = grounded_cells(CellView(Q_g), H, b)
        mask = Function(Q_g, name="dhdt_mask")
        mask.dat.data[:] = cov.dat.data_ro * grounded.dat.data_ro
        area = assemble(mask * dx(mesh))
        if area <= 0.0:
            raise SystemExit("run_eigendec.py --likelihood map: no grounded cell Smith covers")
        L = float(meta.get("dhdt_filter_l", 0.0))
        test = TestFunction(Q)

        def filter_form(out_, src):
            return (out_ * test + Constant(L ** 2) * inner(fd.grad(out_), fd.grad(test))
                    - src * test) * dx

        mask_s = None
        if L > 0.0:
            mask_s = Function(Q, name="dhdt_mask_filtered")
            fd.solve(filter_form(mask_s, mask) == 0, mask_s, solver_parameters=FILTER_SP)
            mask_s.dat.data[:] = np.maximum(mask_s.dat.data_ro, 1e-3)
        out.dhdt = SimpleNamespace(
            weight=w, sigma=float(os.environ.get("ISMIP7_DHDT_SIGMA", "0.1")),
            dt=float(os.environ.get("ISMIP7_DHDT_DT", "1.0")), mask=mask, area=area,
            filter_l=L, filter_form=filter_form, mask_s=mask_s,
            net_sigma=float(meta.get("dhdt_net_sigma", 0.0)),
            total_area=assemble(Constant(1.0) * dx(mesh)))
        record["terms"].append("dhdt")
        record.update({"dhdt_weight": w, "dhdt_sigma": out.dhdt.sigma, "dhdt_dt": out.dhdt.dt,
                       "dhdt_filter_l": L})
        if out.dhdt.net_sigma > 0.0:
            record["terms"].append("dhdt_net")
            record["dhdt_net_sigma"] = out.dhdt.net_sigma
    return out


def dhdt_observations(ctx, terms, u, fc_params):
    r"""The dH/dt term's observed variables at the velocity ``u``: the stepped
    thickness, and the filtered masked tendency when the term is filtered.
    Each comes from an EquationSolver, so a recording manager tapes them and
    the Gauss-Newton action differentiates through them; with the manager
    stopped they are plain solves on new Functions. Empty without the term."""
    dh = terms.dhdt
    if dh is None:
        return []
    H = ctx["h"]
    Q_g = H.function_space()
    dt = Constant(dh.dt)
    h_next = Function(Q_g, name="h_next")
    EquationSolver(transport_form(h_next, H, TestFunction(Q_g), u, dt, Constant(0.0)) == 0,
                   h_next, solver_parameters=TRANSPORT_SP,
                   form_compiler_parameters=fc_params).solve()
    if dh.mask_s is None:
        return [h_next]
    r_s = Function(dh.mask_s.function_space(), name="dhdt_resid_filtered")
    EquationSolver(dh.filter_form(r_s, dh.mask * (h_next - H) / dt) == 0, r_s,
                   solver_parameters=FILTER_SP, form_compiler_parameters=fc_params).solve()
    return [r_s, h_next]


def dhdt_R_inv(terms, misfit_scale, tau):
    r"""The dH/dt term's second derivative with respect to its observed
    variables (``dhdt_observations``), applied to the directions ``tau``:
    the filtered tendency's weight ``misfit_scale w / area * mask /
    (m_s sigma)^2``, or the pointwise one on the stepped thickness, and the
    net term, rank one in the stepped thickness through the masked
    integral of its tendency."""
    dh = terms.dhdt
    out = []
    c = misfit_scale * dh.weight / dh.area / dh.sigma ** 2
    if dh.mask_s is not None:
        t_r, t_h = tau
        Q = t_r.function_space()
        out.append(assemble(Constant(c) * dh.mask / dh.mask_s ** 2 * t_r * TestFunction(Q) * dx,
                            tensor=fd.Cofunction(Q.dual())))
    else:
        t_h, = tau
    Q_g = t_h.function_space()
    psi = TestFunction(Q_g)
    form = Constant(0.0) * t_h * psi * dx
    if dh.mask_s is None:
        form = Constant(c / dh.dt ** 2) * dh.mask * t_h * psi * dx
    R_h = assemble(form, tensor=fd.Cofunction(Q_g.dual()))
    if dh.net_sigma > 0.0:
        # 0.5 misfit_scale (m c)^2 with m = int(mask (h - H)) / (dt A) and
        # c = A rho / sigma_net: rank one in the stepped thickness
        c_net = dh.total_area * RHO_GT / dh.net_sigma
        k = misfit_scale * c_net ** 2 / (dh.dt * dh.total_area) ** 2
        R_h.dat.data[:] += (k * float(assemble(dh.mask * t_h * dx))
                            * assemble(dh.mask * psi * dx).dat.data_ro)
    out.append(R_h)
    return out


def likelihood_note(record, meta):
    r"""What the spreads computed under the likelihood ``record`` bound, for
    the MAP whose recorded settings are ``meta``."""
    if record["likelihood"] == "map":
        return ("posterior from the MAP's own objective (terms: "
                + ", ".join(record["terms"]) + "); a truncated spectrum bounds it from above")
    note = "posterior from the velocity term alone"
    if (float(meta.get("dhdt_weight", 0.0)) > 0.0
            or float(meta.get("dhdt_reweight", 0.0)) > 0.0):
        note += ("; the inversion also scored dH/dt (dhdt_weight > 0 or dhdt_reweight > 0), "
                 "so these spreads bound its posterior from above")
    return note


def build_problem(qoi="discharge", likelihood="map"):
    r"""The posterior's ingredients at the MAP ``ISMIP7_INVERSION`` names: the
    forward loaded as a projection loads it, the inversion's data term, one
    taped solve at the recorded state, the Gauss-Newton action on it, the
    prior's covariance and precision, and the quantities (``qoi``:
    ``discharge`` or ``none``) with their values, gradients and prior
    variances. ``likelihood`` ``map`` is the MAP's own objective
    (``map_objective_terms``), ``velocity`` its velocity chi-square alone.
    Collective; the eigensolve and the direct posterior solve both start
    from it."""
    comm = MPI.COMM_WORLD
    # Every Hessian action must leave memory as it found it, or a solve of a
    # few hundred actions ends out of memory hours in (2592549): refuse a
    # tlm_adjoint that keeps the actions' tangent-linear copies.
    attached, _ = tangent_linear_release_check(comm)
    if max(attached) > 0:
        raise SystemExit(
            "the tlm_adjoint this job imports keeps each Gauss-Newton action's tangent-linear "
            f"copies ({attached} still attached after three actions of a small problem), so "
            "the solve would run out of memory; install the pin (antarctica/scripts/"
            "batch_runners/install_deps.sh --check-tlm-adjoint) or put it first on PYTHONPATH")

    ctx = simulation.setup_model()
    mesh, z, sparams = ctx["mesh"], ctx["z"], ctx["sparams"]
    phi = ctx["phi"]
    alpha = ctx["alpha"]
    if alpha is None:
        raise SystemExit("run_eigendec.py: the MAP must carry the sqrt friction control "
                         "(alpha); log/exp controls are not wired yet")
    # The objective and prior the MAP records (setup_model keeps only the
    # attributes a forward needs; the eigendecomposition needs the rest).
    with fd.CheckpointFile(os.environ["ISMIP7_INVERSION"], "r") as chk:
        meta = {k: chk.get_attr("/", k) for k in
                ("misfit_scale", "log_vel_weight", "log_vel_eps", "prior_form",
                 "prior_sigma_alpha", "prior_rho_theta", "prior_sigma_phi", "prior_rho",
                 "dhdt_weight", "dhdt_filter_l", "dhdt_net_sigma", "dhdt_reweight",
                 "dhdt_reweight_scale", "dhdt_speed_ratio")
                if chk.has_attr("/", k)}
    if str(meta.get("prior_form", "")) != "bilaplacian":
        raise SystemExit(f"run_eigendec.py expects a bilaplacian-prior MAP, got "
                         f"prior_form={meta.get('prior_form')!r}")
    build_F = ctx["build_F"]
    if build_F is None:
        raise SystemExit("run_eigendec.py needs a residual friction law (budd)")
    fc_params = {"quadrature_degree": 4}
    Q = phi.function_space()
    V = VectorFunctionSpace(mesh, "CG", 1)

    # The data term, as the inversion built it: chi^2 on the MEaSUREs errors,
    # observed nodes only, times the MAP's recorded misfit scale, plus the
    # log-speed term at its recorded weight. Residual against u_MAP: zero.
    vel_fn = simulation.find_file(os.path.join(obs_data_root(), "velocity"), "*.nc")
    err = icepack.interpolate((rasterio.open(f"netcdf:{vel_fn}:ERRX"),
                               rasterio.open(f"netcdf:{vel_fn}:ERRY")), V, fillvalue=0.0)
    emag = Function(Q).interpolate(sqrt(err[0] ** 2 + err[1] ** 2))
    observed = np.asarray(emag.dat.data_ro > 0.0)
    obs_mask = Function(Q, name="obs_mask").assign(1.0)
    obs_mask.dat.data[:] = observed.astype(float)
    floor = float(os.environ.get("ISMIP7_SIGMA_U_FLOOR", "3.0"))
    sig_ux = Function(Q).interpolate(max_value(abs(err[0]), Constant(floor)))
    sig_uy = Function(Q).interpolate(max_value(abs(err[1]), Constant(floor)))
    sig_ux.dat.data[~observed] = 1e4
    sig_uy.dat.data[~observed] = 1e4
    area_val = assemble(obs_mask * dx(mesh))
    misfit_scale = float(meta.get("misfit_scale", 1.0))
    log_w = float(meta.get("log_vel_weight", 0.0))
    log_eps = Constant(float(meta.get("log_vel_eps", 1.0)))
    u_map = Function(V, name="u_map").assign(z.subfunctions[0])
    PETSc.Sys.Print(f"Likelihood: misfit_scale {misfit_scale:g}, log-speed weight {log_w:g}, "
                    f"sigma floor {floor:g} m/yr, {comm.allreduce(int(observed.sum()))} observed nodes")

    # R^-1: the data term's second derivative with respect to the velocity,
    # as a bilinear form on the mixed state. The term is
    #   misfit_scale * 0.5 / area * mask * sum_i (u_i - d_i)^2 / sigma_i^2
    # plus, with a log-speed weight, log_w * 0.5 / area * mask * r^2 with
    # r = ln((|u| + eps) / (|d| + eps)), whose Gauss-Newton part has the
    # gradient of r at the MAP velocity in it.
    Zs = z.function_space()
    v_u = fd.split(TestFunction(Zs))[0]
    sp_map = sqrt(u_map[0] ** 2 + u_map[1] ** 2 + Constant(1e-12))
    dr_map = u_map / (sp_map * (sp_map + log_eps))

    terms = map_objective_terms(ctx, meta, Q, likelihood)
    vel_weight = terms.vel_weight if terms.vel_weight is not None else Constant(1.0)
    PETSc.Sys.Print(f"Likelihood ({likelihood}): " + ", ".join(
        f"{k_} {v_}" for k_, v_ in terms.record.items() if k_ != "likelihood"))

    def R_inv_action(*tau):
        t_u = fd.split(tau[0])[0]
        form = (Constant(misfit_scale) / area_val * obs_mask * vel_weight
                * (t_u[0] * v_u[0] / sig_ux ** 2 + t_u[1] * v_u[1] / sig_uy ** 2) * dx(mesh))
        if log_w > 0.0:
            form = form + (Constant(log_w) / area_val * obs_mask
                           * inner(dr_map, t_u) * inner(dr_map, v_u) * dx(mesh))
        out = [assemble(form, tensor=fd.Cofunction(Zs.dual()))]
        if terms.dhdt is not None:
            out += dhdt_R_inv(terms, misfit_scale, tau[1:])
        return out

    # Record the forward at the MAP: one taped solve from the converged state.
    adjoint_sparams = {k: v for k, v in sparams.items() if k != "snes_atol"}
    # z is converged already: the taped solve must exit at iteration 0 on an
    # absolute tolerance (setup_model set that on its own SNES, not in the
    # dict) with the step-size exit live, or it grinds 200 iterations at the
    # residual floor and reports divergence.
    # build_F() with no argument is the residual setup_model solved, with the
    # controls in it as the forward holds them: alpha^2 as the friction and
    # phi through the MAP's own map (on floating ice only under a floating
    # fluidity control). Passing phi itself skips that map, so the recorded
    # equation is one the MAP state does not satisfy and its derivative with
    # respect to phi is not the inversion's.
    F = build_F()
    with assemble(F, form_compiler_parameters=fc_params).dat.vec_ro as _rv:
        fnorm0 = float(_rv.norm())
    tape_sparams = dict(sparams)
    tape_sparams.update({"snes_atol": 100.0 * max(fnorm0, 1e-300),
                         "snes_stol": 1e-8, "snes_max_it": 50})
    PETSc.Sys.Print(f"Taped solve from the converged state: ||F||={fnorm0:.3e}, "
                    f"snes_atol={tape_sparams['snes_atol']:.3e}")
    reset_manager()
    start_manager()
    F = build_F()
    EquationSolver(F == 0, z, solver_parameters=tape_sparams,
                   adjoint_solver_parameters=adjoint_sparams,
                   form_compiler_parameters=fc_params).solve()
    u_s = fd.split(z)[0]
    X = [z] + dhdt_observations(ctx, terms, u_s, fc_params)
    qois = []
    if qoi == "discharge":
        for name, form in discharge_quantities(ctx, u_s):
            Jq = Functional(name=name)
            Jq.assign(form)
            qois.append((name, Jq))
    stop_manager()
    H = CachedGaussNewton(X, R_inv_action)
    controls = [alpha, phi]

    # The prior precision A M^-1 A per control from the MAP's own record.
    test, trial = TestFunction(Q), TrialFunction(Q)
    coeffs = {
        "alpha": bilaplacian_coeffs(float(meta["prior_sigma_alpha"]), float(meta["prior_rho_theta"])),
        "phi": bilaplacian_coeffs(float(meta["prior_sigma_phi"]), float(meta["prior_rho"])),
    }
    fac = {"ksp_type": "preonly", "pc_type": "cholesky", "pc_factor_mat_solver_type": "mumps"}
    A_mat = {k: assemble(prior_operator_form(trial, test, *c)) for k, c in coeffs.items()}
    A_solver = {k: fd.LinearSolver(m, solver_parameters=fac) for k, m in A_mat.items()}
    M_solver = fd.LinearSolver(assemble(inner(trial, test) * dx), solver_parameters=fac)
    PETSc.Sys.Print("Prior (bi-Laplacian A M^-1 A): " + ", ".join(
        f"{k}: delta {c[0]:.3e} gamma {c[1]:.3e}" for k, c in coeffs.items()))

    def covariance(g_alpha, g_phi):
        """Gamma g = A^-1 M A^-1 g for an assembled (dual) pair."""
        out = []
        for k, g in (("alpha", g_alpha), ("phi", g_phi)):
            x = Function(Q)
            A_solver[k].solve(x, g)
            y = Function(Q)
            A_solver[k].solve(y, assemble(inner(x, test) * dx))
            out.append(y)
        return out

    def precision(v_alpha, v_phi):
        """Gamma^-1 v = A M^-1 A v, returned as an assembled dual pair."""
        out = []
        for k, v in (("alpha", v_alpha), ("phi", v_phi)):
            Av = assemble(prior_operator_form(v, test, *coeffs[k]))
            x = Function(Q)
            M_solver.solve(x, Av)
            out.append(assemble(prior_operator_form(x, test, *coeffs[k])))
        return out

    # The quantities' gradients and prior variances, before the (long)
    # eigensolve, so they are on disk even if it does not finish.
    with fd.CheckpointFile(os.environ["ISMIP7_INVERSION"], "r") as chk:
        dhdt_weight = float(chk.get_attr("/", "dhdt_weight")) if chk.has_attr("/", "dhdt_weight") else 0.0
    qoi_rows = []
    grads = []
    if qois:
        values = [float(Jq) for _, Jq in qois]
        dJs = compute_gradient([Jq for _, Jq in qois], controls)
        for (name, _), value, (ga, gp) in zip(qois, values, dJs):
            ya, yp = covariance(ga, gp)
            prior_var = dual_dot(ga, ya) + dual_dot(gp, yp)
            grads.append((ga, gp, ya, yp))
            qoi_rows.append({"name": name, "value_gt_yr": value,
                             "prior_sigma_gt_yr": float(np.sqrt(max(prior_var, 0.0))),
                             "prior_variance": prior_var})
            PETSc.Sys.Print(f"  {name:28s} {value:9.1f} Gt/yr, prior sigma "
                            f"{np.sqrt(max(prior_var, 0.0)):8.2f} Gt/yr")
    return SimpleNamespace(
        comm=comm, ctx=ctx, mesh=mesh, z=z, sparams=sparams, alpha=alpha, phi=phi, Q=Q, H=H,
        controls=controls, R_inv_action=R_inv_action, covariance=covariance, precision=precision,
        qois=qois, grads=grads, qoi_rows=qoi_rows, misfit_scale=misfit_scale, log_w=log_w,
        floor=floor, coeffs=coeffs, meta=meta, dhdt_weight=dhdt_weight, lc=ctx["lc"],
        out_dir=simulation.RESULTS_DIR, likelihood=terms.record,
        observe=lambda zc: [zc] + dhdt_observations(ctx, terms, fd.split(zc)[0], fc_params))


def make_server(problem, label="the eigensolve"):
    r"""The operator service: rank 0 names an operation on a global control
    vector (``alpha`` then ``phi``), every rank applies it. ``H_dual`` is the
    Gauss-Newton action on a primal vector, ``Ginv`` the prior precision on a
    primal vector, ``G`` the prior covariance on a dual one; each returns a
    global vector, ``stop`` ends the other ranks' service loop. Returns the
    service with ``pack``/``unpack`` between function pairs and global
    vectors, the vector lengths ``n`` (one control) and ``N`` (both), the
    Hessian action counter and the clock the action log reads."""
    comm, Q, H, controls = problem.comm, problem.Q, problem.H, problem.controls
    n = len(to_global(problem.phi))
    N = 2 * n
    va, vp = Function(Q), Function(Q)
    n_actions = [0]
    t_start = perf_counter()

    def unpack(x):
        from_global(x[:n], va)
        from_global(x[n:], vp)

    def pack(fa, fp):
        return np.concatenate([to_global(fa), to_global(fp)])

    def serve(op, x=None):
        """Collective: rank 0 names the operation and the vector; all apply."""
        op, x = comm.bcast((op, x), root=0)
        if op == "stop":
            return None
        unpack(x)
        if op == "H_dual":
            t0 = perf_counter()
            ddJ = H.action(controls, [va, vp])
            out = pack(ddJ[0], ddJ[1])
            n_actions[0] += 1
            rss = rss_mb()
            PETSc.Sys.Print(f"  Hessian action {n_actions[0]}: {perf_counter() - t0:.1f} s "
                            f"({(perf_counter() - t_start) / 3600:.2f} h into {label}); "
                            f"resident {comm.allreduce(rss, op=MPI.SUM) / 1024:.1f} GB over "
                            f"{comm.size} ranks, {comm.allreduce(rss, op=MPI.MAX) / 1024:.2f} GB max")
        elif op == "Ginv":
            out = pack(*problem.precision(va, vp))
        elif op == "G":
            # the covariance acts on duals (assembled gradients)
            ga, gp = fd.Cofunction(Q.dual()), fd.Cofunction(Q.dual())
            from_global(x[:n], ga)
            from_global(x[n:], gp)
            out = pack(*problem.covariance(ga, gp))
        return out

    return SimpleNamespace(serve=serve, pack=pack, unpack=unpack, n=n, N=N, n_actions=n_actions)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modes", type=int, default=int(os.environ.get("ISMIP7_EIG_MODES", "40")))
    ap.add_argument("--tag", default=os.environ.get("ISMIP7_RUN_TAG", ""))
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="ARPACK's relative accuracy on the eigenvalues (default 1e-6; 0 is machine precision)")
    ap.add_argument("--maxiter", type=int, default=None,
                    help="ARPACK restarts allowed; the solve then takes at most 1 + ncv + maxiter (ncv - modes) "
                         "Hessian actions and keeps the modes converged by then, so a job can be sized to its "
                         "wall-clock limit (default: unbounded)")
    ap.add_argument("--ncv", type=int, default=None,
                    help="Lanczos vectors per restart cycle (default 2 modes + 1, at least 20)")
    ap.add_argument("--qoi", choices=("discharge", "none"), default="discharge",
                    help="quantities whose posterior spread is computed (default: discharge)")
    ap.add_argument("--likelihood", choices=("map", "velocity"), default="map",
                    help="the data term the Hessian is of: map, the MAP's own objective (the "
                         "velocity chi-square with its recorded dH/dt re-weighting, the dH/dt "
                         "term and the net term), or velocity, the velocity chi-square alone")
    ap.add_argument("--check", action="store_true",
                    help="check the derivatives and stop: the total discharge's adjoint gradient "
                         "against a central difference, and the Gauss-Newton action d.Hd against "
                         "the R^-1 norm of the velocity change a central difference gives along d")
    args = ap.parse_args()
    if args.maxiter is not None and args.maxiter < 1:
        ap.error(f"--maxiter must be at least 1, got {args.maxiter}")
    P = build_problem(args.qoi, args.likelihood)
    if args.check:
        if not P.qois:
            raise SystemExit("--check needs the discharge quantities (--qoi discharge)")
        # The Gauss-Newton action first: the checks then move the controls,
        # and the cached Gauss-Newton refuses to act once they have moved.
        check_gauss_newton(P.ctx, P.z, P.sparams, P.alpha, P.phi, P.H, P.controls, P.grads[0],
                           P.R_inv_action, P.observe)
        check_gradient(P.ctx, P.z, P.sparams, P.alpha, P.phi, P.qois[0], P.grads[0])
        PETSc.Sys.Print("Derivative checks done; run without --check for the eigendecomposition.")
        return
    comm, qois, grads, qoi_rows = P.comm, P.qois, P.grads, P.qoi_rows
    mesh, Q, dhdt_weight, misfit_scale, log_w, floor, coeffs, meta = (
        P.mesh, P.Q, P.dhdt_weight, P.misfit_scale, P.log_w, P.floor, P.coeffs, P.meta)
    S = make_server(P)
    serve, n, N, n_actions = S.serve, S.n, S.N, S.n_actions
    k = min(args.modes, N - 2)
    ncv = min(args.ncv or arpack_ncv(k), N)
    if comm.rank == 0:
        budget = (f"at most {arpack_action_bound(k, ncv, args.maxiter)} Hessian actions "
                  f"({args.maxiter} restarts)" if args.maxiter is not None else "restarts unbounded")
        PETSc.Sys.Print(f"ARPACK: {k} leading generalised eigenpairs of H_GN v = lambda Gamma^-1 v "
                        f"({N} control dofs); {ncv} Lanczos vectors, tol {args.tol:g}, {budget}")
        try:
            vals, vecs, converged = leading_eigenpairs(
                lambda x: serve("H_dual", np.asarray(x, dtype=float)),
                lambda x: serve("Ginv", np.asarray(x, dtype=float)),
                lambda x: serve("G", np.asarray(x, dtype=float)),
                N, k, tol=args.tol, maxiter=args.maxiter, ncv=ncv)
        except BaseException:
            traceback.print_exc()
            sys.stderr.flush()
            comm.Abort(1)
        serve("stop")
    else:
        while serve(None) is not None:
            pass
        vals, vecs, converged = None, None, None
    vals, converged = comm.bcast((vals, converged), root=0)
    if not converged:
        PETSc.Sys.Print(f"ARPACK reached its {args.maxiter} restarts with {vals.size} of {k} modes converged; "
                        f"the spreads below use those (still upper bounds)")
    k = int(vals.size)
    PETSc.Sys.Print(f"Leading eigenvalues (lambda), {n_actions[0]} Hessian actions: "
                    + " ".join(f"{v:.3e}" for v in vals[:10]))
    PETSc.Sys.Print(f"  modes with lambda > 1 (data-dominated): {int((vals > 1.0).sum())} of {k}; "
                    f"min {vals.min():.3e}")

    lc = P.lc
    sfx = f"_{args.tag}" if args.tag else ""
    out_dir = P.out_dir
    os.makedirs(out_dir, exist_ok=True)
    eig_fn = os.path.join(out_dir, f"eigendec_{lc}{sfx}.h5")
    with fd.CheckpointFile(eig_fn, "w") as chk:
        chk.save_mesh(mesh)
        projections = np.zeros((len(grads), k))
        for i in range(k):
            x = comm.bcast(vecs[:, i] if comm.rank == 0 else None, root=0)
            ma, mp = Function(Q, name=f"mode_{i:04d}_alpha"), Function(Q, name=f"mode_{i:04d}_phi")
            from_global(x[:n], ma)
            from_global(x[n:], mp)
            chk.save_function(ma)
            chk.save_function(mp)
            for q, (ga, gp, _, _) in enumerate(grads):
                projections[q, i] = dual_dot(ga, ma) + dual_dot(gp, mp)
        for (name, _), (_, _, ya, yp) in zip(qois[:2], grads[:2]):
            chk.save_function(ya, name=f"sensitivity_{name}_alpha")
            chk.save_function(yp, name=f"sensitivity_{name}_phi")
        chk.set_attr("/", "n_modes", k)
        chk.set_attr("/", "map", os.path.basename(os.environ["ISMIP7_INVERSION"]))
    if comm.rank == 0:
        np.savetxt(os.path.join(out_dir, f"eigenvalues_{lc}{sfx}.txt"), vals)
        with open(os.path.join(out_dir, f"eigendec_{lc}{sfx}.json"), "w") as f:
            json.dump({"map": os.environ["ISMIP7_INVERSION"], "modes": k,
                       "modes_requested": args.modes, "converged": bool(converged),
                       "arpack": {"tol": args.tol, "ncv": ncv, "maxiter": args.maxiter,
                                  "hessian_actions": n_actions[0]},
                       "misfit_scale": misfit_scale, "log_vel_weight": log_w,
                       "sigma_floor": floor, "likelihood": P.likelihood,
                       "prior": {k_: list(c) for k_, c in coeffs.items()},
                       "prior_record": {k_: meta.get(k_) for k_ in
                                        ("prior_sigma_alpha", "prior_rho_theta", "prior_sigma_phi", "prior_rho")},
                       "eigenvalues": [float(v) for v in vals]}, f, indent=2)
    PETSc.Sys.Print(f"Saved: {eig_fn}")

    if qoi_rows:
        lam = np.clip(vals, 0.0, None)
        PETSc.Sys.Print(f"Discharge uncertainty from the inferred friction and fluidity "
                        f"({k} modes; likelihood terms {', '.join(P.likelihood['terms'])}):")
        for q, row in enumerate(qoi_rows):
            post_var = posterior_variance(row["prior_variance"], projections[q], lam)
            share = lam / (1.0 + lam) * projections[q] ** 2
            row.update({"posterior_sigma_gt_yr": float(np.sqrt(max(post_var, 0.0))),
                        "posterior_variance": post_var,
                        "variance_reduction_by_mode": [float(v) for v in share],
                        "last_mode_share_of_prior": float(share[-1] / row["prior_variance"])
                        if row["prior_variance"] > 0 else 0.0})
            PETSc.Sys.Print(f"  {row['name']:28s} {row['value_gt_yr']:9.1f} Gt/yr  sigma prior "
                            f"{row['prior_sigma_gt_yr']:8.2f} -> posterior "
                            f"{row['posterior_sigma_gt_yr']:8.2f} Gt/yr")
        if comm.rank == 0:
            with open(os.path.join(out_dir, f"uncertainty_{lc}{sfx}.json"), "w") as f:
                json.dump({"map": os.environ["ISMIP7_INVERSION"], "modes": k,
                           "modes_requested": args.modes, "converged": bool(converged),
                           "likelihood": P.likelihood,
                           "map_dhdt_weight": dhdt_weight,
                           "note": likelihood_note(P.likelihood, meta),
                           "eigenvalues": [float(v) for v in vals],
                           "quantities": qoi_rows}, f, indent=2)


def _solve_at(ctx, z, sparams, alpha, phi, a0, p0, da, dp, eps):
    """Re-solve the momentum balance, untaped, at the controls (a0, p0) + eps (da, dp),
    from the MAP state, on a copy of it: the recorded z is state-locked by the
    Gauss-Newton object."""
    params = dict(sparams)
    params.update({"snes_rtol": 1e-10, "snes_atol": 0.0, "snes_max_it": 100})
    alpha.assign(a0 + eps * da)
    phi.assign(p0 + eps * dp)
    zc = z.copy(deepcopy=True)
    fd.solve(ctx["build_F"](z_c=zc) == 0, zc, solver_parameters=params,
             form_compiler_parameters={"quadrature_degree": 4})
    return zc


def check_gauss_newton(ctx, z, sparams, alpha, phi, H, controls, grad, R_inv_action, observe):
    """d . H_GN d against the R^-1 norm of the change in every observed
    variable (``observe``: the state, and the dH/dt term's stepped thickness
    and filtered tendency under the MAP's likelihood) along d, from a central
    difference of two untaped solves, d = Gamma g of the total."""
    _, _, ya, yp = grad
    comm = alpha.function_space().mesh().comm
    scale = max(float(np.max(np.abs(ya.dat.data_ro), initial=0.0)),
                float(np.max(np.abs(yp.dat.data_ro), initial=0.0)))
    scale = comm.allreduce(scale, op=MPI.MAX)
    if scale == 0.0:
        PETSc.Sys.Print("Gauss-Newton check skipped: zero direction")
        return
    da = Function(ya.function_space()).assign(ya / scale)
    dp = Function(yp.function_space()).assign(yp / scale)
    from time import perf_counter
    t0 = perf_counter()
    Hd = H.action(controls, [da, dp])
    PETSc.Sys.Print(f"  one Gauss-Newton action: {perf_counter() - t0:.1f} s "
                    f"(the eigensolve needs about three per mode)")
    dHd = dual_dot(Hd[0], da) + dual_dot(Hd[1], dp)
    a0, p0 = alpha.copy(deepcopy=True), phi.copy(deepcopy=True)
    eps = 1e-3
    z_plus = _solve_at(ctx, z, sparams, alpha, phi, a0, p0, da, dp, +eps)
    z_minus = _solve_at(ctx, z, sparams, alpha, phi, a0, p0, da, dp, -eps)
    alpha.assign(a0)
    phi.assign(p0)
    X_plus, X_minus = observe(z_plus), observe(z_minus)
    tau = [Function(xp.function_space()).assign((xp - xm) / (2.0 * eps))
           for xp, xm in zip(X_plus, X_minus)]
    norm_fd = sum(dual_dot(r, t) for r, t in zip(R_inv_action(*tau), tau))
    PETSc.Sys.Print(f"Gauss-Newton check: d.Hd {dHd:.6e}, |J d|^2_R from a central difference "
                    f"{norm_fd:.6e}, ratio {dHd / norm_fd if norm_fd else float('nan'):.6f}")


def check_gradient(ctx, z, sparams, alpha, phi, qoi, grad):
    """Central difference of the total discharge along Gamma g against the
    adjoint pairing g . Gamma g, from two untaped solves at alpha +- eps d."""
    name, _ = qoi
    ga, gp, ya, yp = grad
    a0, p0 = alpha.copy(deepcopy=True), phi.copy(deepcopy=True)
    scale = max(float(np.max(np.abs(ya.dat.data_ro), initial=0.0)),
                float(np.max(np.abs(yp.dat.data_ro), initial=0.0)))
    scale = alpha.function_space().mesh().comm.allreduce(scale, op=MPI.MAX)
    eps = 1e-3 / scale if scale > 0 else 0.0
    values = []
    for sign in (+1.0, -1.0):
        zc = _solve_at(ctx, z, sparams, alpha, phi, a0, p0, ya, yp, sign * eps)
        values.append(float(assemble(dict(discharge_quantities(ctx, fd.split(zc)[0]))[name])))
    alpha.assign(a0)
    phi.assign(p0)
    fd_slope = (values[0] - values[1]) / (2.0 * eps)
    adj_slope = dual_dot(ga, ya) + dual_dot(gp, yp)
    PETSc.Sys.Print(f"Gradient check ({name}): central difference {fd_slope:.6e}, adjoint "
                    f"{adj_slope:.6e}, ratio {fd_slope / adj_slope if adj_slope else float('nan'):.6f}")


if __name__ == "__main__":
    main()
