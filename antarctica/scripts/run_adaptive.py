r"""Run a forward experiment with Úa-style mesh adaptation between segments.

    python antarctica/scripts/run_adaptive.py \
        --driver antarctica/scripts/control/run.py --experiment-name ctrl2015_cesm2_waccm_adapt \
        --launcher "mpiexec -n 12" --t-start 2015 --t-end 2115 --adapt-every 10 \
        --initial-iterations 2 [--tag adapt]

Úa adapts the mesh before the first run-step (`AdaptMeshInitial`, iterated up
to `AdaptMeshMaxIterations`) and then every `AdaptMeshTimeInterval`. This
orchestrator does the same with the existing pieces, unchanged:

  segment k : the ordinary driver runs from its restart to t_k
              (ISMIP7_RESTART / ISMIP7_T_END / ISMIP7_RUN_TAG), writing
              <experiment>_<lc>_final.h5 as it always does
  adapt     : antarctica/scripts/adapt_mesh.py moves that checkpoint onto a
              new mesh (size field from the ISMIP7_ADAPT_* configuration)
  repeat.

The initial adaptation runs the driver for ZERO years first (t_end = t_start)
so the diagnostic velocity exists on the starting mesh for the strain-rate
criteria, adapts with --rebuild-aref (the apparent-mass-balance reference is
rebuilt exactly on the new mesh), and iterates until the element count
changes by less than --until-change (Úa's
AdaptMeshUntilChangeInNumberOfElementsLessThan) or --initial-iterations is
reached. ISMIP7_APPARENT_MB must be set in the environment for a balanced
control, as for any run.
"""
import argparse
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RESULTS = os.path.join(HERE, "..", "results")


def sh(cmd, env):
    print(f"\n$ {cmd}\n", flush=True)
    rc = subprocess.call(cmd, shell=True, env=env, cwd=REPO)
    if rc != 0:
        raise SystemExit(f"command failed ({rc}): {cmd}")


def checkpoint_year(chk):
    r"""The model year the driver recorded in ``chk``, or None if it has none.

    Firedrake writes ``t_yr`` as a plain HDF5 root attribute, so this reads it
    without importing Firedrake and keeps this orchestrator dependency-free.
    """
    import h5py
    with h5py.File(chk, "r") as f:
        t = f["/"].attrs.get("t_yr")
    return None if t is None else float(t)


def n_cells(msh):
    r"""Element count of a gmsh 2.2 file, from its $Elements block."""
    n = 0
    with open(msh) as f:
        for line in f:
            if line.startswith("$Elements"):
                total = int(next(f))
                for _ in range(total):
                    parts = next(f).split()
                    if parts[1] == "2":
                        n += 1
                break
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--driver", required=True)
    ap.add_argument("--experiment-name", required=True,
                    help="the driver's experiment_name (with --tag applied), e.g. ctrl2015_cesm2_waccm_adapt")
    ap.add_argument("--launcher", default="mpiexec -n 4")
    ap.add_argument("--adapt-launcher", default="mpiexec -n 1",
                    help="launcher for the adapt step only. Single-rank by default because "
                         "ISMIP7_ADAPT_TRANSFER=project (the recommended transfer) refuses to "
                         "run on more than one rank, and the remesh is serial gmsh either way. "
                         "Raise it to --launcher's rank count under the shipped default "
                         "transfer (interpolate), which has no such restriction")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--t-start", type=float, required=True)
    ap.add_argument("--t-end", type=float, required=True)
    ap.add_argument("--adapt-every", type=float, required=True, help="years (Úa AdaptMeshTimeInterval)")
    ap.add_argument("--initial-iterations", type=int, default=1, help="Úa AdaptMeshMaxIterations (0 = no initial adaptation)")
    ap.add_argument("--until-change", type=int, default=0, help="Úa AdaptMeshUntilChangeInNumberOfElementsLessThan")
    ap.add_argument("--tag", default=None, help="forwarded as ISMIP7_RUN_TAG")
    ap.add_argument("--restart", default=os.environ.get("ISMIP7_RESTART"), help="start from this checkpoint instead of a cold start")
    ap.add_argument("--lc", default=os.environ.get("ISMIP7_LC", "2500"))
    args = ap.parse_args()

    env = dict(os.environ)
    if args.tag:
        env["ISMIP7_RUN_TAG"] = args.tag
    env.pop("ISMIP7_AUTO_RESUME", None)          # segments are explicit here
    final = os.path.join(RESULTS, f"{args.experiment_name}_{args.lc}_final.h5")
    restart = args.restart
    k = 0

    def run_to(t_end, restart):
        e = dict(env)
        e["ISMIP7_T_END"] = f"{t_end:g}"
        if restart:
            e["ISMIP7_RESTART"] = restart
        else:
            e.pop("ISMIP7_RESTART", None)
        sh(f"{args.launcher} {shlex.quote(args.python)} -u {shlex.quote(args.driver)}", e)
        if not os.path.exists(final):
            raise SystemExit(f"driver did not write {final}")
        # `final` exists from the previous segment too, so its presence proves
        # nothing. The driver exits 0 after a stall (simulation.py saves the
        # last converged year and returns), and adapting a stalled state then
        # re-attempting the same years loops silently to the end of the
        # timeline. Compare the year it actually recorded against the one asked
        # for; the segment counter is not evidence.
        t_got = checkpoint_year(final)
        if t_got is None:
            print(f"WARNING: {final} has no t_yr attribute; cannot confirm the "
                  f"segment reached {t_end:g}", flush=True)
        elif t_got < t_end - 1e-6:
            raise SystemExit(
                f"segment stopped at t_yr={t_got:g}, short of {t_end:g}: the "
                f"driver saved its last converged year and exited. Adapting "
                f"this state would re-attempt the same years. See {final}."
            )
        return final

    def adapt(chk, k, rebuild):
        out = chk.replace("_final.h5", f"_adapt{k}.h5")
        flag = " --rebuild-aref" if rebuild else ""
        sh(f"{args.adapt_launcher} {shlex.quote(args.python)} -u {shlex.quote(os.path.join(HERE, 'adapt_mesh.py'))} "
           f"{shlex.quote(chk)} --out-checkpoint {shlex.quote(out)}{flag}", env)
        return out

    t0 = time.time()
    # --- initial adaptation: Úa AdaptMeshInitial, iterated --------------
    if args.initial_iterations > 0:
        chk = run_to(args.t_start, restart)          # zero years: diagnostics on the start mesh
        prev = None
        for it in range(args.initial_iterations):
            k += 1
            restart = adapt(chk, k, rebuild=True)
            chk = run_to(args.t_start, restart)      # zero years again on the new mesh
            import glob
            newest = max(glob.glob(os.path.join(HERE, "..", "mesh", "*_adapt*.msh")), key=os.path.getmtime)
            n = n_cells(newest)
            print(f"\n=== initial adaptation {it + 1}/{args.initial_iterations}: {n} cells ===", flush=True)
            if prev is not None and abs(n - prev) < args.until_change:
                print("=== element count settled; stopping the initial iteration ===", flush=True)
                break
            prev = n
        restart = chk
    # --- the run, adapting every --adapt-every years --------------------
    t = args.t_start
    while t < args.t_end - 1e-9:
        t_next = min(t + args.adapt_every, args.t_end)
        chk = run_to(t_next, restart)
        t = t_next
        if t < args.t_end - 1e-9:
            k += 1
            restart = adapt(chk, k, rebuild=False)
    print(f"\n=== done: {args.t_start:g} -> {args.t_end:g} with {k} adaptations in {(time.time() - t0) / 3600:.2f} h ===", flush=True)


if __name__ == "__main__":
    main()
