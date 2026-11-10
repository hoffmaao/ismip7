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
The likelihood is the velocity term, as in the eigensolve.
"""
import argparse
import json
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
from icepack2_tools.laplace import load_modes, posterior_variance_cg  # noqa: E402

STATE_KEYS = ("x", "r", "p", "alphas", "betas", "rho", "seed", "projections")


def sigma(variance):
    return float(np.sqrt(max(float(variance), 0.0)))


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
    ap.add_argument("--wall-hours", type=float, default=None,
                    help="stop starting steps this long after the script began")
    args = ap.parse_args()
    if args.iterations < 1:
        ap.error(f"--iterations must be at least 1, got {args.iterations}")
    t_job = perf_counter()

    P = run_eigendec.build_problem("discharge")
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
        meta = None
        if root:
            with open(eig_fn[:-3] + ".json") as f:
                meta = json.load(f)
        meta = comm.bcast(meta, root=0)
        if os.path.basename(str(meta["map"])) != os.path.basename(os.environ["ISMIP7_INVERSION"]):
            raise SystemExit(f"{eig_fn} holds the modes of {meta['map']}, "
                             f"this MAP is {os.environ['ISMIP7_INVERSION']}")
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

    out_fn = os.path.join(P.out_dir, f"posterior_{P.lc}{sfx}.json")
    results = {}
    if root and args.resume and os.path.exists(out_fn):
        with open(out_fn) as f:
            results = json.load(f).get("quantities", {})

    def write(results):
        with open(out_fn, "w") as f:
            json.dump({"map": os.environ["ISMIP7_INVERSION"], "likelihood": "velocity",
                       "map_dhdt_weight": P.dhdt_weight,
                       "note": ("bracket of the Laplace posterior variance from the velocity "
                                "likelihood: lower = A-norm of the conjugate-gradient iterate, "
                                "upper = Gauss-Radau quadrature with its node at 1; the zero-step "
                                "upper bound is the truncated-spectrum variance of the seed modes"),
                       "modes": int(lam.size) if lam is not None else 0,
                       "eigendec": eig_fn if lam is not None else None,
                       "lower_eigenvalue": 1.0,
                       "written": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                       "quantities": results}, f, indent=2)

    def over_budget():
        return args.wall_hours is not None and perf_counter() - t_job > 3600.0 * args.wall_hours

    for name in wanted:
        q = names.index(name)
        ga, gp, _, _ = P.grads[q]
        row = P.qoi_rows[q]
        g = S.pack(ga, gp)
        state_fn = os.path.join(P.out_dir, f"posterior_cg_{P.lc}{sfx}_{name}.npz")
        summary = None
        if root:
            state = None
            if args.resume and os.path.exists(state_fn):
                with np.load(state_fn) as z:
                    state = {k: z[k] for k in STATE_KEYS}
                PETSc.Sys.Print(f"{name}: continuing from {len(state['alphas'])} steps in {state_fn}")
            record = dict(results.get(name, {}))
            record.update({"value_gt_yr": row["value_gt_yr"], "prior_sigma_gt_yr": row["prior_sigma_gt_yr"],
                           "prior_variance": row["prior_variance"]})
            t_q = perf_counter()

            def callback(step):
                record.update({"iterations": step["iteration"],
                               "sigma_lower_gt_yr": sigma(step["lower"]),
                               "sigma_upper_gt_yr": sigma(step["upper"]),
                               "variance_lower": step["lower"], "variance_upper": step["upper"],
                               "rho": step["rho"], "seconds": perf_counter() - t_q})
                record.setdefault("lower_by_step", []).append(step["lower"])
                record.setdefault("upper_by_step", []).append(step["upper"])
                record.setdefault("rho_by_step", []).append(step["rho"])
                results[name] = record
                write(results)
                PETSc.Sys.Print(f"  {name} step {step['iteration']}: sigma {sigma(step['lower']):.1f} to "
                                f"{sigma(step['upper']):.1f} Gt/yr, rho {step['rho']:.3e}")

            try:
                out = posterior_variance_cg(
                    lambda x: S.serve("H_dual", np.asarray(x, dtype=float)),
                    lambda x: S.serve("Ginv", np.asarray(x, dtype=float)),
                    lambda x: S.serve("G", np.asarray(x, dtype=float)),
                    g, args.iterations, modes=V, eigenvalues=lam, state=state, callback=callback,
                    stop=lambda k: over_budget())
            except BaseException:
                traceback.print_exc()
                sys.stderr.flush()
                comm.Abort(1)
            S.serve("stop")
            np.savez(state_fn, **{k: np.asarray(out[k]) for k in STATE_KEYS})
            if state is None:
                record["lower_by_step"] = [float(v) for v in out["lower"]]
                record["upper_by_step"] = [float(v) for v in out["upper"]]
                record["rho_by_step"] = [float(v) for v in out["rho"]]
            record.update({"iterations": out["iterations"],
                           "hessian_actions": int(record.get("hessian_actions", 0)) + out["actions"],
                           "truncation_sigma_gt_yr": sigma(out["seed"] + out["rho"][0]),
                           "seed_variance": out["seed"],
                           "variance_lower": float(out["lower"][-1]), "variance_upper": float(out["upper"][-1]),
                           "sigma_lower_gt_yr": sigma(out["lower"][-1]),
                           "sigma_upper_gt_yr": sigma(out["upper"][-1]),
                           "rho": float(out["rho"][-1]), "seconds": perf_counter() - t_q,
                           "state": state_fn})
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
