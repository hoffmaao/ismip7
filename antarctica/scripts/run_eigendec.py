#!/usr/bin/env python3
r"""Leading eigenmodes of the prior-preconditioned Gauss-Newton Hessian of a MAP.

The Laplace approximation of Recinos et al. (2023) and fenics_ice: at the MAP
the posterior covariance is ``(H_GN + Gamma^-1)^-1`` with ``H_GN`` the
Gauss-Newton Hessian of the data term and ``Gamma`` the prior covariance.
Its departure from the prior lives in the generalised eigenpairs

    H_GN v_i = lambda_i Gamma^-1 v_i,

ordered by lambda_i: a mode with lambda_i >> 1 is constrained by the data,
one with lambda_i << 1 is known no better than the prior. Those pairs are what
this script computes and saves, and what a forward's QoI sensitivity is
projected onto (``sigma_Q^2 = sigma_prior^2 - sum_i (g.v_i)^2 lambda_i/(1+lambda_i)``).

The problem is the MAP's own: ``simulation.setup_model()`` loads the MAP
(``ISMIP7_INVERSION``) exactly as a forward does -- its mesh, geometry, friction
control (sqrt, exp or log), sub-element scheme and the mixed state it was
accepted at -- and the data term is the inversion's (the per-datum chi^2 on
the MEaSUREs errors with the MAP's recorded ``misfit_scale`` and log-speed
weight), evaluated against the MAP's own velocity so the residual is zero and
the Hessian is Gauss-Newton: positive semi-definite by construction. The
prior is the MAP's recorded bi-Laplacian (sigma, rho per control).

The tape holds one solve at the physical exponents from the converged state
(the direct forward of 1 Oct 2026); the n=1->3 ladder inside the tape is what
made the earlier Hessian indefinite. Run it under the full mixed-Jacobian
MUMPS the inversion differentiates through:

    ISMIP7_INVERSION=<map.h5> ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=full_mumps \
        mpiexec -n 4 python antarctica/scripts/run_eigendec.py [--modes 40]

Output beside the MAP's results directory: ``eigendec_<lc>[_<tag>].h5`` (the
modes as alpha/phi pairs, Gamma^-1-orthonormal), ``eigenvalues_<lc>[_<tag>].txt``
and ``eigendec_<lc>[_<tag>].json`` (the prior and likelihood the modes belong
to). The eigensolve is ARPACK on rank 0 over the global control vector; every
rank serves the collective Hessian and prior actions.
"""
import argparse
import json
import os
import sys

import numpy as np
import firedrake as fd
from firedrake import (Constant, Function, TestFunction,
                       TrialFunction, assemble, dx, inner, ln, max_value, sqrt,
                       VectorFunctionSpace)
from firedrake.petsc import PETSc
from mpi4py import MPI
from scipy.sparse.linalg import LinearOperator, eigsh
from tlm_adjoint.firedrake import (CachedHessian, EquationSolver, Functional,
                                   reset_manager, start_manager, stop_manager)

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

import icepack  # noqa: E402
import rasterio  # noqa: E402
import simulation  # noqa: E402
from icepack2_tools.prior import (bilaplacian_coeffs,  # noqa: E402
                                  prior_operator_form)
from icepack2_tools.runconfig import obs_data_root  # noqa: E402


def to_global(f):
    with f.dat.vec_ro as v:
        scatter, x_seq = PETSc.Scatter.toAll(v)
        scatter.scatter(v, x_seq, mode=PETSc.Scatter.Mode.FORWARD)
        out = x_seq.array.copy()
        scatter.destroy(); x_seq.destroy()
    return out


def from_global(arr, f):
    with f.dat.vec_wo as v:
        x_seq = PETSc.Vec().createSeq(len(arr), comm=PETSc.COMM_SELF)
        x_seq.array[:] = arr
        scatter, _ = PETSc.Scatter.toAll(v)
        scatter.scatter(x_seq, v, mode=PETSc.Scatter.Mode.REVERSE)
        scatter.destroy(); x_seq.destroy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modes", type=int, default=int(os.environ.get("ISMIP7_EIG_MODES", "40")))
    ap.add_argument("--tag", default=os.environ.get("ISMIP7_RUN_TAG", ""))
    args = ap.parse_args()
    comm = MPI.COMM_WORLD

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
                 "prior_sigma_alpha", "prior_rho_theta", "prior_sigma_phi", "prior_rho")
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

    def data_term(u):
        chi2 = (0.5 / area_val * obs_mask
                * ((u[0] - u_map[0]) ** 2 / sig_ux ** 2 + (u[1] - u_map[1]) ** 2 / sig_uy ** 2))
        J = Constant(misfit_scale) * chi2 * dx(mesh)
        if log_w > 0.0:
            sp = sqrt(u[0] ** 2 + u[1] ** 2 + Constant(1e-12))
            sp0 = sqrt(u_map[0] ** 2 + u_map[1] ** 2 + Constant(1e-12))
            J = J + Constant(log_w) * (0.5 / area_val * obs_mask
                                       * ln((sp + log_eps) / (sp0 + log_eps)) ** 2) * dx(mesh)
        return J

    # Record the forward at the MAP: one taped solve from the converged state.
    adjoint_sparams = {k: v for k, v in sparams.items() if k != "snes_atol"}
    # z is converged already: the taped solve must exit at iteration 0 on an
    # absolute tolerance (setup_model set that on its own SNES, not in the
    # dict) with the step-size exit live, or it grinds 200 iterations at the
    # residual floor and reports divergence.
    F = build_F(phi_c=phi)
    with assemble(F, form_compiler_parameters=fc_params).dat.vec_ro as _rv:
        fnorm0 = float(_rv.norm())
    tape_sparams = dict(sparams)
    tape_sparams.update({"snes_atol": 100.0 * max(fnorm0, 1e-300),
                         "snes_stol": 1e-8, "snes_max_it": 50})
    PETSc.Sys.Print(f"Taped solve from the converged state: ||F||={fnorm0:.3e}, "
                    f"snes_atol={tape_sparams['snes_atol']:.3e}")
    reset_manager()
    start_manager()
    F = build_F(phi_c=phi)
    EquationSolver(F == 0, z, solver_parameters=tape_sparams,
                   adjoint_solver_parameters=adjoint_sparams,
                   form_compiler_parameters=fc_params).solve()
    u_s = fd.split(z)[0]
    J = Functional(name="J")
    J.assign(data_term(u_s))
    stop_manager()
    PETSc.Sys.Print(f"J at the MAP: {float(J):.3e} (zero residual by construction)")
    H = CachedHessian(J)
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
            x = Function(Q); A_solver[k].solve(x, g)
            y = Function(Q); A_solver[k].solve(y, assemble(inner(x, test) * dx))
            out.append(y)
        return out

    def precision(v_alpha, v_phi):
        """Gamma^-1 v = A M^-1 A v, returned as an assembled dual pair."""
        out = []
        for k, v in (("alpha", v_alpha), ("phi", v_phi)):
            Av = assemble(prior_operator_form(v, test, *coeffs[k]))
            x = Function(Q); M_solver.solve(x, Av)
            out.append(assemble(prior_operator_form(x, test, *coeffs[k])))
        return out

    n = len(to_global(phi))
    N = 2 * n
    va, vp = Function(Q), Function(Q)
    n_actions = [0]

    def unpack(x):
        from_global(x[:n], va); from_global(x[n:], vp)

    def pack(fa, fp):
        return np.concatenate([to_global(fa), to_global(fp)])

    def serve(op, x=None):
        """Collective: rank 0 names the operation and the vector; all apply."""
        op, x = comm.bcast((op, x), root=0)
        if op == "stop":
            return None
        unpack(x)
        if op == "H_dual":
            from time import perf_counter
            t0 = perf_counter()
            _, _, ddJ = H.action(controls, [va, vp])
            out = pack(ddJ[0], ddJ[1])
            n_actions[0] += 1
            if n_actions[0] % 10 == 1:
                PETSc.Sys.Print(f"  Hessian action {n_actions[0]}: {perf_counter() - t0:.1f} s")
        elif op == "Ginv":
            out = pack(*precision(va, vp))
        elif op == "G":
            # the covariance acts on duals (assembled gradients)
            ga, gp = fd.Cofunction(Q.dual()), fd.Cofunction(Q.dual())
            from_global(x[:n], ga); from_global(x[n:], gp)
            out = pack(*covariance(ga, gp))
        return out

    k = min(args.modes, N - 2)
    if comm.rank == 0:
        PETSc.Sys.Print(f"ARPACK: {k} leading generalised eigenpairs of H_GN v = lambda Gamma^-1 v "
                        f"({N} control dofs)")
        H_op = LinearOperator((N, N), matvec=lambda x: serve("H_dual", np.asarray(x, dtype=float)), dtype=float)
        M_op = LinearOperator((N, N), matvec=lambda x: serve("Ginv", np.asarray(x, dtype=float)), dtype=float)
        Minv_op = LinearOperator((N, N), matvec=lambda x: serve("G", np.asarray(x, dtype=float)), dtype=float)
        vals, vecs = eigsh(H_op, k=k, M=M_op, Minv=Minv_op, which="LA")
        order = np.argsort(-vals)
        vals, vecs = vals[order], vecs[:, order]
        serve("stop")
    else:
        while serve(None) is not None:
            pass
        vals, vecs = None, None
    vals = comm.bcast(vals, root=0)
    PETSc.Sys.Print("Leading eigenvalues (lambda): " + " ".join(f"{v:.3e}" for v in vals[:10]))
    PETSc.Sys.Print(f"  modes with lambda > 1 (data-dominated): {int((vals > 1.0).sum())} of {k}; "
                    f"min {vals.min():.3e}")

    lc = ctx["lc"]
    sfx = f"_{args.tag}" if args.tag else ""
    out_dir = simulation.RESULTS_DIR
    os.makedirs(out_dir, exist_ok=True)
    eig_fn = os.path.join(out_dir, f"eigendec_{lc}{sfx}.h5")
    with fd.CheckpointFile(eig_fn, "w") as chk:
        chk.save_mesh(mesh)
        for i in range(k):
            x = comm.bcast(vecs[:, i] if comm.rank == 0 else None, root=0)
            ma, mp = Function(Q, name=f"mode_{i:04d}_alpha"), Function(Q, name=f"mode_{i:04d}_phi")
            from_global(x[:n], ma); from_global(x[n:], mp)
            chk.save_function(ma); chk.save_function(mp)
        chk.set_attr("/", "n_modes", k)
        chk.set_attr("/", "map", os.path.basename(os.environ["ISMIP7_INVERSION"]))
    if comm.rank == 0:
        np.savetxt(os.path.join(out_dir, f"eigenvalues_{lc}{sfx}.txt"), vals)
        with open(os.path.join(out_dir, f"eigendec_{lc}{sfx}.json"), "w") as f:
            json.dump({"map": os.environ["ISMIP7_INVERSION"], "modes": k,
                       "misfit_scale": misfit_scale, "log_vel_weight": log_w,
                       "sigma_floor": floor, "prior": {k_: list(c) for k_, c in coeffs.items()},
                       "prior_record": {k_: meta.get(k_) for k_ in
                                        ("prior_sigma_alpha", "prior_rho_theta", "prior_sigma_phi", "prior_rho")},
                       "eigenvalues": [float(v) for v in vals]}, f, indent=2)
    PETSc.Sys.Print(f"Saved: {eig_fn}")


if __name__ == "__main__":
    main()
