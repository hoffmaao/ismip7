#!/usr/bin/env python3
r"""Posterior spread of the grounded discharge by direct solves at the MAP.

The eigendecomposition (``run_eigendec.py``) bounds a quantity's Laplace
posterior variance from above by its leading modes. At 2 km the velocity
data dominate the prior in far more directions than an eigensolve can hold
(the hundredth eigenvalue of the production MAP is still 3.5e5), so the
truncated spread keeps falling with every mode and never settles. This script
brackets the variance of each quantity instead: conjugate gradients on
``(H + Gamma^-1) x = g`` preconditioned with the prior covariance
(``icepack2_tools.laplace.posterior_variance_cg``), where ``g . x_k`` rises to
the variance from below and Gauss-Radau quadrature on the same recurrence
bounds it from above. The saved modes seed the solve with the truncated
posterior and are deflated from it, so the bracket starts at the
eigensolve's own upper bound and tightens with every Hessian action, one per
step, at the cost the eigensolve logged per action.

    ISMIP7_INVERSION=<map.h5> mpiexec -n 32 python antarctica/scripts/run_posterior_variance.py \
        --tag <tag> --quantities discharge_total,grounding_line_flux_total --iterations 200

The problem, the quantities, their gradients and the operator service are
``run_eigendec.build_problem`` and ``run_eigendec.make_server``; the modes
come from ``eigendec_<lc>_<tag>.h5`` beside the results (``--eigendec``,
``--modes`` for fewer; 0 solves from the prior). Each quantity's bracket is
written to ``posterior_<lc>_<tag>.json`` after every step, with the vectors
to continue from in ``posterior_cg_<lc>_<tag>_<quantity>.npz``; ``--resume``
continues those, and ``--wall-hours`` stops the job cleanly before its limit.
The likelihood (``--likelihood``, the MAP's own objective by default) must be
the one the seed modes were computed under; the record beside them is checked.
A resume must run under the likelihood its saved vectors and bracket were
computed under, and state saved with no record counts as the velocity term.
"""
import argparse
import json
import math
import os
import sys
import traceback
from datetime import datetime, timezone
from time import perf_counter

import numpy as np
from firedrake.petsc import PETSc

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

import run_eigendec  # noqa: E402
from icepack2_tools.laplace import (load_modes, posterior_variance_cg,  # noqa: E402
                                    same_likelihood)

STATE_KEYS = ("x", "r", "p", "alphas", "betas", "rho", "seed", "projections")
LOWER_EIGENVALUE = 1.0
# the record of output written before the likelihood was recorded, which ran the velocity term alone
VELOCITY_ONLY = {"likelihood": "velocity", "terms": ["velocity"]}


def sigma(variance):
    return float(np.sqrt(max(float(variance), 0.0)))


def or_abort(comm, read):
    r"""``read()`` on this rank; any failure prints its traceback and aborts
    every rank, since the others wait in a collective."""
    try:
        return read()
    except BaseException:
        traceback.print_exc()
        sys.stderr.flush()
        comm.Abort(1)


def read_json(fn):
    with open(fn) as f:
        return json.load(f)


def replace_atomically(fn, write):
    r"""``write(f)`` into a sibling temporary file, then moved onto ``fn``, so a
    kill mid-write leaves the previous complete file."""
    tmp = fn + ".tmp"
    with open(tmp, "wb") as f:
        write(f)
    os.replace(tmp, fn)


def check_operator(meta, P, eig_fn):
    r"""Stop unless the eigensolve that wrote ``meta`` ran the live likelihood
    and prior, so its modes are eigenvectors of this operator."""
    live = {"sigma_floor": P.floor, "misfit_scale": P.misfit_scale, "log_vel_weight": P.log_w}
    live.update({f"prior {k}": list(c) for k, c in P.coeffs.items()})
    saved = {k: meta.get(k) for k in ("sigma_floor", "misfit_scale", "log_vel_weight")}
    saved.update({f"prior {k}": (meta.get("prior") or {}).get(k) for k in P.coeffs})
    have_l = meta.get("likelihood") or VELOCITY_ONLY
    if not same_likelihood(have_l, P.likelihood):
        raise SystemExit(f"{eig_fn} was computed under the likelihood {have_l}, "
                         f"this run has {P.likelihood}")
    for key, want in live.items():
        have = saved[key]
        want_v = np.atleast_1d(np.asarray(want, dtype=float))
        have_v = None if have is None else np.atleast_1d(np.asarray(have, dtype=float))
        if have_v is None or have_v.shape != want_v.shape or not all(
                math.isclose(h, w, rel_tol=1e-9, abs_tol=0.0) for h, w in zip(have_v, want_v)):
            raise SystemExit(f"{eig_fn} was computed with {key} {have}, this run has {want}")


def check_resume(z, n_modes, eig_name, lam, likelihood, state_fn):
    r"""The reason the saved state ``z`` cannot continue under the live modes
    and the live likelihood record, or None when it can."""
    have_l = json.loads(str(z["likelihood"])) if "likelihood" in z.files else VELOCITY_ONLY
    if not same_likelihood(have_l, likelihood):
        return f"{state_fn} was built under the likelihood {have_l}, this run has {likelihood}"
    saved = {k: z[k] if k in z.files else None
             for k in ("n_modes", "eigendec", "eigenvalues", "lower_eigenvalue")}
    if saved["n_modes"] is None or int(saved["n_modes"]) != n_modes \
            or str(saved["eigendec"]) != eig_name \
            or float(saved["lower_eigenvalue"]) != LOWER_EIGENVALUE \
            or not np.array_equal(np.asarray(saved["eigenvalues"], dtype=float), lam):
        return (f"{state_fn} was built from {saved['n_modes']} modes of {saved['eigendec']} "
                f"(lower eigenvalue {saved['lower_eigenvalue']}), this run seeds with {n_modes} "
                f"modes of {eig_name!r} (lower eigenvalue {LOWER_EIGENVALUE})")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default=os.environ.get("ISMIP7_RUN_TAG", ""))
    ap.add_argument("--quantities", default="discharge_total,grounding_line_flux_total",
                    help="comma-separated names from the eigensolve's uncertainty file, or all")
    ap.add_argument("--iterations", type=int, default=200,
                    help="conjugate-gradient steps per quantity, one Hessian action each")
    ap.add_argument("--modes", type=int, default=None,
                    help="saved modes to seed and deflate with (default all; 0 for none)")
    ap.add_argument("--eigendec", default=None,
                    help="the eigendecomposition file (default eigendec_<lc>_<tag>.h5 beside the results)")
    ap.add_argument("--resume", action="store_true",
                    help="continue each quantity from its saved vectors")
    ap.add_argument("--likelihood", choices=("map", "velocity"), default="map",
                    help="the data term (run_eigendec.py --likelihood); the seed modes must "
                         "have been computed under the same one")
    ap.add_argument("--wall-hours", type=float, default=None,
                    help="stop starting steps this long after the script began")
    args = ap.parse_args()
    if args.iterations < 1:
        ap.error(f"--iterations must be at least 1, got {args.iterations}")
    t_job = perf_counter()

    P = run_eigendec.build_problem("discharge", args.likelihood)
    comm = P.comm
    root = comm.rank == 0
    S = run_eigendec.make_server(P, label="the posterior solves")
    sfx = f"_{args.tag}" if args.tag else ""
    names = [name for name, _ in P.qois]
    wanted = names if args.quantities == "all" else args.quantities.split(",")
    unknown = [w for w in wanted if w not in names]
    if unknown:
        raise SystemExit(f"run_posterior_variance.py: unknown quantities {unknown}; "
                         f"the MAP offers {names}")

    eig_fn = args.eigendec or os.path.join(P.out_dir, f"eigendec_{P.lc}{sfx}.h5")
    V = lam = None
    if args.modes != 0 and os.path.exists(eig_fn):
        meta_fn = eig_fn[:-3] + ".json"
        if not os.path.exists(meta_fn):
            raise SystemExit(f"{eig_fn} has no record {meta_fn}")
        meta = None
        if root:
            meta = or_abort(comm, lambda: read_json(meta_fn))
        meta = comm.bcast(meta, root=0)
        if os.path.basename(str(meta["map"])) != os.path.basename(os.environ["ISMIP7_INVERSION"]):
            raise SystemExit(f"{eig_fn} holds the modes of {meta['map']}, "
                             f"this MAP is {os.environ['ISMIP7_INVERSION']}")
        check_operator(meta, P, eig_fn)
        modes, copied = load_modes(eig_fn, P.Q, args.modes)
        lam = np.asarray(meta["eigenvalues"][:len(modes)], dtype=float)
        packed = [S.pack(ma, mp) for ma, mp in modes]
        del modes
        if root:
            V = np.array(packed)
        del packed
        PETSc.Sys.Print(f"Seed: {lam.size} modes from {eig_fn} ({'copied' if copied else 'interpolated'}), "
                        f"lambda {lam[0]:.3e} to {lam[-1]:.3e}; the zero-step upper bound is the "
                        f"truncated variance")
    elif args.modes != 0:
        PETSc.Sys.Print(f"No eigendecomposition at {eig_fn}: solving from the prior")
    else:
        PETSc.Sys.Print("Solving from the prior (--modes 0)")

    n_modes = int(lam.size) if lam is not None else 0
    eig_name = os.path.basename(eig_fn) if lam is not None else ""
    lam_saved = lam if lam is not None else np.zeros(0)
    out_fn = os.path.join(P.out_dir, f"posterior_{P.lc}{sfx}.json")
    results = {}
    refusal = None
    if root and args.resume and os.path.exists(out_fn):
        doc = or_abort(comm, lambda: read_json(out_fn))
        have_l = doc.get("likelihood") or VELOCITY_ONLY
        if same_likelihood(have_l, P.likelihood):
            results = doc.get("quantities", {})
        else:
            refusal = f"{out_fn} was written under the likelihood {have_l}, this run has {P.likelihood}"
    refusal = comm.bcast(refusal, root=0)
    if refusal is not None:
        raise SystemExit(refusal)

    def write(results):
        doc = {"map": os.environ["ISMIP7_INVERSION"], "likelihood": P.likelihood,
               "map_dhdt_weight": P.dhdt_weight,
               "note": ("bracket of the Laplace posterior variance under the likelihood terms "
                        + ", ".join(P.likelihood["terms"]) + ": lower = A-norm of the "
                        "conjugate-gradient iterate, upper = Gauss-Radau quadrature with its "
                        "node at 1; the zero-step upper bound is the truncated-spectrum "
                        "variance of the seed modes"),
               "modes": n_modes,
               "eigendec": eig_fn if lam is not None else None,
               "lower_eigenvalue": LOWER_EIGENVALUE,
               "written": datetime.now(timezone.utc).isoformat(timespec="seconds"),
               "quantities": results}
        replace_atomically(out_fn, lambda f: f.write(json.dumps(doc, indent=2).encode()))

    def over_budget():
        return args.wall_hours is not None and perf_counter() - t_job > 3600.0 * args.wall_hours

    for name in wanted:
        q = names.index(name)
        ga, gp, _, _ = P.grads[q]
        row = P.qoi_rows[q]
        g = S.pack(ga, gp)
        state_fn = os.path.join(P.out_dir, f"posterior_cg_{P.lc}{sfx}_{name}.npz")
        summary = state = refusal = None
        if root and args.resume and os.path.exists(state_fn):
            def read_state():
                with np.load(state_fn) as z:
                    refusal = check_resume(z, n_modes, eig_name, lam_saved, P.likelihood, state_fn)
                    return refusal, (None if refusal else {k: z[k] for k in STATE_KEYS})

            refusal, state = or_abort(comm, read_state)
        refusal = comm.bcast(refusal, root=0)
        if refusal is not None:
            raise SystemExit(refusal)
        if root:
            if state is not None:
                PETSc.Sys.Print(f"{name}: continuing from {len(state['alphas'])} steps in {state_fn}")
            record = dict(results.get(name, {}))
            record.update({"value_gt_yr": row["value_gt_yr"], "prior_sigma_gt_yr": row["prior_sigma_gt_yr"],
                           "prior_variance": row["prior_variance"]})
            t_q = perf_counter()

            actions_before = int(record.get("hessian_actions", 0)) if state is not None else 0

            def callback(step):
                replace_atomically(state_fn, lambda f: np.savez(
                    f, n_modes=n_modes, eigendec=eig_name, eigenvalues=lam_saved,
                    lower_eigenvalue=LOWER_EIGENVALUE, likelihood=json.dumps(P.likelihood), **{k: np.asarray(step[k]) for k in STATE_KEYS}))
                rho = float(step["rho"][-1])
                record.update({"iterations": step["iteration"],
                               "hessian_actions": actions_before + step["actions"],
                               "truncation_sigma_gt_yr": sigma(step["seed"] + step["rho"][0]),
                               "seed_variance": step["seed"],
                               "sigma_lower_gt_yr": sigma(step["lower"]),
                               "sigma_upper_gt_yr": sigma(step["upper"]),
                               "variance_lower": float(step["lower"]),
                               "variance_upper": float(step["upper"]),
                               "rho": rho, "seconds": perf_counter() - t_q,
                               "lower_by_step": [float(v) for v in step["lower_by_step"]],
                               "upper_by_step": [float(v) for v in step["upper_by_step"]],
                               "rho_by_step": [float(v) for v in step["rho"]],
                               "state": state_fn})
                results[name] = record
                write(results)
                PETSc.Sys.Print(f"  {name} step {step['iteration']}: sigma {sigma(step['lower']):.1f} to "
                                f"{sigma(step['upper']):.1f} Gt/yr, rho {rho:.3e}")

            try:
                out = posterior_variance_cg(
                    lambda x: S.serve("H_dual", np.asarray(x, dtype=float)),
                    lambda x: S.serve("Ginv", np.asarray(x, dtype=float)),
                    lambda x: S.serve("G", np.asarray(x, dtype=float)),
                    g, args.iterations, modes=V, eigenvalues=lam, lower_eigenvalue=LOWER_EIGENVALUE,
                    state=state, callback=callback,
                    stop=lambda k: over_budget())
            except BaseException:
                traceback.print_exc()
                sys.stderr.flush()
                comm.Abort(1)
            S.serve("stop")
            record.update({"iterations": out["iterations"],
                           "hessian_actions": actions_before + out["actions"],
                           "truncation_sigma_gt_yr": sigma(out["seed"] + out["rho"][0]),
                           "seed_variance": out["seed"],
                           "variance_lower": float(out["lower"][-1]), "variance_upper": float(out["upper"][-1]),
                           "sigma_lower_gt_yr": sigma(out["lower"][-1]),
                           "sigma_upper_gt_yr": sigma(out["upper"][-1]),
                           "seconds": perf_counter() - t_q, "state": state_fn})
            results[name] = record
            write(results)
            summary = {k: record[k] for k in ("iterations", "sigma_lower_gt_yr", "sigma_upper_gt_yr",
                                              "truncation_sigma_gt_yr", "hessian_actions")}
        else:
            while S.serve(None) is not None:
                pass
        summary = comm.bcast(summary, root=0)
        PETSc.Sys.Print(f"{name}: {row['value_gt_yr']:.1f} Gt/yr, sigma prior {row['prior_sigma_gt_yr']:.1f}, "
                        f"truncated {summary['truncation_sigma_gt_yr']:.1f}, bracket "
                        f"{summary['sigma_lower_gt_yr']:.1f} to {summary['sigma_upper_gt_yr']:.1f} Gt/yr "
                        f"after {summary['iterations']} steps ({summary['hessian_actions']} Hessian actions)")
        if comm.bcast(over_budget(), root=0):
            PETSc.Sys.Print(f"Wall budget of {args.wall_hours} h reached; resume with --resume")
            break
    PETSc.Sys.Print(f"Saved: {out_fn}")


if __name__ == "__main__":
    main()
