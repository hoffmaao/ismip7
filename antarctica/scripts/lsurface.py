#!/usr/bin/env python3
r"""L-surface over the two regularisation weights, one inversion per point.

The inversion has two priors, on log-friction (theta, ISMIP7_GAMMA_THETA) and
log-fluidity (phi, ISMIP7_GAMMA_PHI), so the classical L-curve is a surface:
misfit against the two prior energies over a (gamma_theta, gamma_phi) grid.
Each point is an ordinary `submit.sh inversion` run with its own MAP and
timing JSON; nothing here solves anything.

    # queue a 4x4 grid, two decades either side of 1e5, under the sweep
    # tolerance so a point costs ~a node-day at 2 km
    python scripts/lsurface.py submit --law regularized_coulomb \
        ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_MESH=... ISMIP7_BUFFER_M=0

    # once they finish: table, corners of the fixed-gamma slices, figure
    python scripts/lsurface.py harvest results/lsurface/regularized_coulomb_2000

The corner of each slice is the point of maximum curvature of the log-log
L-curve (Hansen & O'Leary 1993). Two slices through the grid give one corner per weight; the surface
plot shows whether they are consistent.
"""
import argparse
import glob
import json
import os
import subprocess
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ANT = os.path.dirname(_HERE)
SUBMIT = os.path.join(_HERE, "batch_runners", "submit.sh")

# Two decades either side of the production weight, four points per axis.
GAMMAS_DEFAULT = "1e3,1e4,1e5,1e6"
# The sweep tolerance: a point need not be converged to the production
# 1e-10, only far enough down the L-curve to place it.
FTOL_DEFAULT = "1e-4"


def law_tag(law):
    return {"regularized_coulomb": "rc", "budd": "budd"}.get(law, law)


def point_paths(out_dir, gamma_theta, gamma_phi):
    r"""(map, json) for one grid point; the names carry both weights."""
    stem = os.path.join(out_dir, f"gt{gamma_theta:g}_gp{gamma_phi:g}")
    return stem + ".h5", stem + ".json"


def submit_command(law, gamma_theta, gamma_phi, out_dir, ftol, passthrough, dry_run=False):
    r"""The submit.sh argv for one point. Pure, so a test can read it."""
    map_out, json_out = point_paths(out_dir, gamma_theta, gamma_phi)
    name = f"ls_{law_tag(law)}_t{np.log10(gamma_theta):g}_p{np.log10(gamma_phi):g}"
    cmd = [SUBMIT, "inversion", "--name", name]
    if dry_run:
        cmd.append("--dry-run")
    cmd += [
        f"ISMIP7_FRICTION={law}",
        f"ISMIP7_GAMMA_THETA={gamma_theta:g}",
        f"ISMIP7_GAMMA_PHI={gamma_phi:g}",
        f"ISMIP7_FTOL={ftol}",
        f"ISMIP7_MAP_OUT={map_out}",
        f"ISMIP7_INVERSION_TIMING_JSON={json_out}",
        *passthrough,
    ]
    return cmd


def cmd_submit(args):
    gt = [float(g) for g in args.gammas.split(",")]
    gp = [float(g) for g in (args.gammas_phi or args.gammas).split(",")]
    lc = next((kv.split("=", 1)[1] for kv in args.passthrough if kv.startswith("ISMIP7_LC=")),
              os.environ.get("ISMIP7_LC", ""))
    out_dir = args.out_dir or os.path.join(_ANT, "results", "lsurface", f"{args.law}_{lc}")
    os.makedirs(out_dir, exist_ok=True)
    for g_t in gt:
        for g_p in gp:
            cmd = submit_command(args.law, g_t, g_p, out_dir, args.ftol, args.passthrough, args.dry_run)
            print(" ".join(cmd) if args.dry_run else f"gamma_theta={g_t:g} gamma_phi={g_p:g}", flush=True)
            subprocess.run(cmd, check=True)
    print(f"{len(gt) * len(gp)} points -> {out_dir}")


# ── harvest ──────────────────────────────────────────────────────────────

def load_points(out_dir):
    r"""One row per finished JSON: weights, and the last evaluation's terms."""
    rows = []
    for fn in sorted(glob.glob(os.path.join(out_dir, "gt*_gp*.json"))):
        with open(fn) as f:
            rec = json.load(f)
        if not rec.get("evaluations"):
            continue
        last = rec["evaluations"][-1]
        rows.append({
            "gamma_theta": float(rec["knobs"]["gamma_theta"]),
            "gamma_phi": float(rec["knobs"]["gamma_phi"]),
            "misfit": float(last["misfit"]),
            "reg_theta": float(last["reg_theta"]),
            "reg_phi": float(last["reg_phi"]),
            "nit": int(rec.get("nit", last.get("eval", 0))),
            "finished": rec.get("phase") == "finished",
            "message": rec.get("message", ""),
        })
    return sorted(rows, key=lambda r: (r["gamma_theta"], r["gamma_phi"]))


def lcurve_corner(gammas, misfit, reg):
    r"""Index of maximum curvature of the log-log L-curve, parametrised by
    log(gamma). Interior points only: curvature needs both neighbours.

    Returns -1 when there are fewer than three points."""
    g = np.log10(np.asarray(gammas, float))
    x = np.log10(np.asarray(misfit, float))
    y = np.log10(np.asarray(reg, float))
    if len(g) < 3:
        return -1
    order = np.argsort(g)
    g, x, y = g[order], x[order], y[order]
    dx, dy = np.gradient(x, g), np.gradient(y, g)
    ddx, ddy = np.gradient(dx, g), np.gradient(dy, g)
    kappa = (dx * ddy - dy * ddx) / np.power(dx * dx + dy * dy, 1.5)
    kappa = np.abs(kappa)
    kappa[[0, -1]] = -np.inf
    return int(order[np.argmax(kappa)])


def slices(rows):
    r"""For each fixed gamma_phi, the L-curve over gamma_theta (misfit vs
    reg_theta), and the converse. Yields (axis, fixed, gammas, misfit, reg, corner)."""
    for fixed_key, var_key, reg_key, axis in (
        ("gamma_phi", "gamma_theta", "reg_theta", "theta"),
        ("gamma_theta", "gamma_phi", "reg_phi", "phi"),
    ):
        for fixed in sorted({r[fixed_key] for r in rows}):
            sl = sorted((r for r in rows if r[fixed_key] == fixed), key=lambda r: r[var_key])
            gam = [r[var_key] for r in sl]
            mis = [r["misfit"] for r in sl]
            reg = [r[reg_key] for r in sl]
            yield axis, fixed, gam, mis, reg, lcurve_corner(gam, mis, reg)


def cmd_harvest(args):
    rows = load_points(args.out_dir)
    if not rows:
        sys.exit(f"no finished points under {args.out_dir}")
    print(f"{'gamma_theta':>12} {'gamma_phi':>10} {'misfit':>12} {'reg_theta':>11} {'reg_phi':>11} {'nit':>4}  status")
    for r in sorted(rows, key=lambda r: (r["gamma_theta"], r["gamma_phi"])):
        print(f"{r['gamma_theta']:12.3g} {r['gamma_phi']:10.3g} {r['misfit']:12.5e} "
              f"{r['reg_theta']:11.4e} {r['reg_phi']:11.4e} {r['nit']:4d}  "
              f"{'finished' if r['finished'] else 'running'}")
    print("\ncorners (max curvature of each log-log slice):")
    for axis, fixed, gam, mis, reg, k in slices(rows):
        other = "gamma_phi" if axis == "theta" else "gamma_theta"
        where = f"gamma_{axis}={gam[k]:g}" if k >= 0 else "n/a (<3 points)"
        print(f"  {other}={fixed:<8g} -> {where}")
    csv = os.path.join(args.out_dir, "lsurface.csv")
    with open(csv, "w") as f:
        f.write("gamma_theta,gamma_phi,misfit,reg_theta,reg_phi,nit,finished\n")
        for r in rows:
            f.write(f"{r['gamma_theta']:g},{r['gamma_phi']:g},{r['misfit']:.8e},"
                    f"{r['reg_theta']:.8e},{r['reg_phi']:.8e},{r['nit']},{int(r['finished'])}\n")
    print(f"\nwrote {csv}")
    if args.plot:
        plot(rows, args.plot)
        print(f"wrote {args.plot}")


def plot(rows, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, want in zip(axes[:2], ("theta", "phi")):
        for axis, fixed, gam, mis, reg, k in slices(rows):
            if axis != want:
                continue
            other = r"$\gamma_\phi$" if want == "theta" else r"$\gamma_\theta$"
            line, = ax.loglog(mis, reg, "o-", label=f"{other}={fixed:g}")
            if k >= 0:
                ax.plot(mis[k], reg[k], "s", ms=11, mfc="none", color=line.get_color())
        ax.set_xlabel("misfit")
        ax.set_ylabel(rf"prior energy, $\{want}$")
        ax.set_title(rf"L-curves over $\gamma_\{want}$ (squares: max curvature)")
        ax.legend(fontsize=8)
    gt = sorted({r["gamma_theta"] for r in rows})
    gp = sorted({r["gamma_phi"] for r in rows})
    z = np.full((len(gp), len(gt)), np.nan)
    for r in rows:
        z[gp.index(r["gamma_phi"]), gt.index(r["gamma_theta"])] = np.log10(r["misfit"])
    ax = axes[2]
    m = ax.pcolormesh(np.log10(gt), np.log10(gp), z, shading="nearest")
    fig.colorbar(m, ax=ax, label=r"$\log_{10}$ misfit")
    ax.set_xlabel(r"$\log_{10}\gamma_\theta$")
    ax.set_ylabel(r"$\log_{10}\gamma_\phi$")
    ax.set_title("L-surface")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("submit", help="queue one inversion per grid point through submit.sh")
    s.add_argument("--law", default="regularized_coulomb", choices=["regularized_coulomb", "budd"])
    s.add_argument("--gammas", default=GAMMAS_DEFAULT, help="gamma_theta values (comma separated)")
    s.add_argument("--gammas-phi", default=None, help="gamma_phi values; defaults to --gammas")
    s.add_argument("--ftol", default=FTOL_DEFAULT)
    s.add_argument("--out-dir", default=None, help="MAPs and JSONs; default results/lsurface/<law>_<lc>")
    s.add_argument("--dry-run", action="store_true", help="print the sbatch lines, start nothing")
    s.add_argument("passthrough", nargs="*", help="KEY=VALUE for submit.sh (mesh, lc, buffer, ...)")
    s.set_defaults(func=cmd_submit)
    h = sub.add_parser("harvest", help="table, slice corners, CSV and figure from the JSONs")
    h.add_argument("out_dir")
    h.add_argument("--plot", default=None, help="write this PNG")
    h.set_defaults(func=cmd_harvest)
    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
