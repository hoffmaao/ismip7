#!/usr/bin/env python3
r"""Run an ISMIP7 forward with a calving law from hoffmaao/calving.

    mpiexec -n N python antarctica/scripts/forward_calving.py --law vonmises \
        --law-param sigma_max_fl=0.2 [--experiment control] [driver args]

``--law`` is one of ``calving/laws.py``'s registered laws (fixed, velocity,
position, thickness, vonmises, vonmises_strain, hfb, or a ``--law-module``
plugin), found in ``--calving-dir`` (default ``$CALVING_DIR``, a checkout of
github.com/hoffmaao/calving). Its ``rate(model, t)`` is evaluated every
transport step on the forward's live dual state (``simulation.LiveCalvingState``)
and handed to the shared level set as the ``prescribed`` rate, so the front
the forward runs is the front the law was tuned on (``calving/tune_greene.py``
against the Greene et al. 2022 ice-front record).

The experiment driver is wrapped, not copied: its ``run_simulation`` call is
intercepted to put the law object in ``ctx`` and everything else (forcing,
restart, tags, auto-resume, ISMIP7_* knobs) stays the driver's own. Leave
``ISMIP7_CALVING`` unset; the law object owns the front.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

ap = argparse.ArgumentParser(add_help=False)
ap.add_argument("--law", default="vonmises")
ap.add_argument("--law-param", action="append", default=[], metavar="KEY=VALUE")
ap.add_argument("--law-module", default=None)
ap.add_argument("--calving-dir", default=os.environ.get("CALVING_DIR"),
                help="checkout of github.com/hoffmaao/calving")
ap.add_argument("--experiment", default="control",
                choices=("control", "ocx", "ssp126", "ssp370", "ssp585",
                         "hist"))
ap.add_argument("--esm", default="cesm_waccm", choices=("cesm_waccm", "mri_esm2"))
args, rest = ap.parse_known_args()
sys.argv = [sys.argv[0]] + rest
if not args.calving_dir or not os.path.isdir(args.calving_dir):
    raise SystemExit("name the hoffmaao/calving checkout with --calving-dir "
                     "or CALVING_DIR")

sys.path.insert(0, args.calving_dir)
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import laws  # noqa: E402  (hoffmaao/calving)

if args.law_module:
    laws.load_module(args.law_module)
law = laws.make(args.law, 1, **laws.parse_params(args.law_param))

import simulation  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402

_real_run_simulation = simulation.run_simulation


def _run_simulation_with_law(ctx, *a, **kw):
    ctx["calving_law"] = law
    PETSc.Sys.Print(f"  External calving law: {law.describe()} "
                    f"(from {args.calving_dir})")
    return _real_run_simulation(ctx, *a, **kw)


simulation.run_simulation = _run_simulation_with_law

if args.experiment == "control":
    sys.path.insert(0, os.path.join(HERE, "control"))
    import run as driver  # control/run.py
    driver.run_simulation = _run_simulation_with_law
    driver.main()
elif args.experiment == "ocx":
    sys.path.insert(0, os.path.join(HERE, "projections"))
    import ocx as driver
    driver.run_simulation = _run_simulation_with_law
    driver.main()
else:
    import runpy
    import experiment
    experiment.run_simulation = _run_simulation_with_law
    sub = "historical" if args.experiment == "hist" else "projections"
    name = (f"{args.esm}.py" if args.experiment == "hist"
            else f"{args.experiment}_{args.esm}.py")
    runpy.run_path(os.path.join(HERE, sub, name), run_name="__main__")
