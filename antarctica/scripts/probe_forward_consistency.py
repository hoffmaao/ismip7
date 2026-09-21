#!/usr/bin/env python3
r"""How sensitive is the re-solve of a MAP's own velocity to each forward-only
argument of the residual?

`check_budd_map.py --forward` measures the distance between the forward's
cold-start diagnostic and the velocity the inversion saved; a MAP inverted
under the forward's residual returns about 1e-7. Both codes build the residual
through `dual_friction.build_rc_residual` and differ only in what they pass.
This runs the same re-solve under each variant, one fresh process per variant
because `simulation.setup_model()` holds global state:

    as_is        the forward as shipped
    nref_none    N_ref=None, the inversion's own call (ISMIP7_BUDD_NREF)
    no_drag      ocean_drag and the soft speed limiter off
    alpha_inv    composite alpha at the inversion's 1e-2 (Budd forward: 1e-4)
    all          every one of the above together

On a fresh 32 km Budd MAP N_ref=None and alpha_inv were inert and the drags
moved the answer; on the 14 September adaptive-mesh MAPs the drags were the whole
difference (FORWARD_RUN_READINESS.md section 4).

    python antarctica/scripts/probe_forward_consistency.py MAP.h5 \
        ISMIP7_MESH=... ISMIP7_LC=32000 ISMIP7_LC_COARSE=320000 [-n 4]
"""
import argparse
import json
import os
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))

VARIANTS = {
    "as_is": {},
    "nref_none": {"ISMIP7_BUDD_NREF": "none"},
    "no_drag": {"ISMIP7_OCEAN_DRAG": "0", "ISMIP7_U_LIM": "0", "ISMIP7_RC_EPS_TAUC": "0"},
    "alpha_inv": {"ISMIP7_COMPOSITE_ALPHA": "1e-2"},
}
VARIANTS["all"] = {k: v for d in VARIANTS.values() for k, v in d.items()}

_CHILD = r"""
import json, os, sys
sys.path.insert(0, %r)
from check_budd_map import forward_check
rel, umean = forward_check(%r)
print("PROBE_RESULT " + json.dumps({"rel_l2": rel, "mean_speed": umean}))
"""


def run_variant(name, extra_env, map_path, nranks, base_env):
    env = dict(base_env)
    env.update(extra_env)
    cmd = ["mpiexec", "-n", str(nranks)] if nranks > 1 else []
    cmd += [sys.executable, "-c", _CHILD % (_HERE, map_path)]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    hits = [ln for ln in proc.stdout.splitlines() if ln.startswith("PROBE_RESULT ")]
    if proc.returncode != 0 or not hits:
        tail = "\n".join((proc.stderr or proc.stdout).splitlines()[-12:])
        return {"error": f"exit {proc.returncode}", "tail": tail}
    return json.loads(hits[0][len("PROBE_RESULT "):])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("map")
    ap.add_argument("env", nargs="*", help="KEY=VALUE for every variant (mesh, lc, data roots)")
    ap.add_argument("-n", "--nranks", type=int, default=1)
    ap.add_argument("--only", default=None, help="comma-separated variant names")
    ap.add_argument("--json", default=None, help="write the table here as well")
    args = ap.parse_args(argv)

    base = dict(os.environ)
    base.setdefault("OMP_NUM_THREADS", "1")
    for kv in args.env:
        k, v = kv.split("=", 1)
        base[k] = v
    names = args.only.split(",") if args.only else list(VARIANTS)

    results = {}
    print(f"{'variant':<10} {'rel L2':>10} {'mean speed':>11}  overrides")
    for name in names:
        res = run_variant(name, VARIANTS[name], os.path.abspath(args.map), args.nranks, base)
        results[name] = res
        if "error" in res:
            print(f"{name:<10} {'FAILED':>10} {'':>11}  {res['error']}\n{res['tail']}")
        else:
            print(f"{name:<10} {res['rel_l2']:10.3e} {res['mean_speed']:11.2f}  "
                  + " ".join(f"{k}={v}" for k, v in VARIANTS[name].items()))
        sys.stdout.flush()
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"map": os.path.abspath(args.map), "env": args.env, "results": results}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
