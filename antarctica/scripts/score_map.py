#!/usr/bin/env python3
r"""Score a MAP by the flux it carries across its own grounding line.

A velocity misfit says how well an inversion fits the observations it was
given; it does not say whether the resulting state will hold its ice. The
quantity that decides that is the grounding-line discharge, and the honest
way to read it is against the flux the OBSERVED velocity carries across the
same facets with the same thickness:

.. code::

    ratio = Q(u_model) / Q(u_obs)

A ratio of one means the inverted velocity delivers ice to the shelves at the
observed rate, so a control run's grounded budget is closed by the same
discharge the real ice sheet has. Below one, grounded ice thickens; the
2500 m MAP sits at 0.75 and its control gains ~1 mm SLE/yr of volume above
flotation (``antarctica/scripts/region_budget.py``).

The same ratio is reported per band of observed speed, because a
sigma-normalised velocity misfit buys its fit where the observational error
is smallest -- the slow interior -- and starves the intermediate tributaries
that carry most of the discharge. That is the bias ISSM's logarithmic misfit
(``ISMIP7_LOG_VEL_WEIGHT``) exists to remove, and this script is how the two
objectives are compared.

Periodic MAP checkpoints carry the controls but not a velocity, so the
diagnostic solve is re-run here through ``simulation.setup_model``: the same
operator the forward would use, so the flux scored is the flux a forward run
would apply.

    python antarctica/scripts/score_map.py MAP.h5 [MAP2.h5 ...]

The run environment (``ISMIP7_LC``, ``ISMIP7_FRICTION``, ``ISMIP7_N_FLOW``,
``ISMIP7_GEOMETRY_SPACE``, ``ISMIP7_MESH``) must match the inversion's, as it
must for any forward that loads the MAP.
"""

import argparse
import os
import sys

import firedrake as fd
from firedrake import Constant, Function

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(_ROOT))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from region_budget import crossing_flux, regions

BANDS = [0.0, 100.0, 500.0, 1500.0, 1e9]
BAND_LABELS = ["< 100", "100 - 500", "500 - 1500", "> 1500"]
#: Rignot et al. 2019 grounding-line discharge.
OBSERVED_DISCHARGE = (2050.0, 100.0)


def score(path):
    from simulation import setup_model
    os.environ["ISMIP7_INVERSION"] = os.path.abspath(path)
    ctx = setup_model()
    mesh, h, b = ctx["mesh"], ctx["h"], ctx["b"]
    u = ctx["z"].subfunctions[0]
    u_obs = ctx["u_obs"]
    state = {"mesh": mesh, "thickness": h, "bed": b, "velocity": u}
    reg = regions(state)
    Q0 = reg["ice"].function_space()
    open_water = Function(Q0)
    open_water.dat.data[:] = 1.0 - reg["ice"].dat.data_ro
    reg["open"] = open_water

    def discharge(velocity, speed=None, lo=None, hi=None):
        """Grounding line plus grounded front: every facet grounded ice leaves
        through, optionally restricted to one observed-speed band."""
        return (crossing_flux(state, reg, "grounded", "floating",
                              speed, lo, hi, velocity)
                + crossing_flux(state, reg, "grounded", "open",
                                speed, lo, hi, velocity))

    q_m = discharge(u)
    q_o = discharge(u_obs)
    sp = Function(Q0).interpolate(
        fd.sqrt(fd.dot(u_obs, u_obs) + Constant(1e-12)))
    rows = []
    for lab, lo, hi in zip(BAND_LABELS, BANDS[:-1], BANDS[1:]):
        m = discharge(u, sp, lo, hi)
        o = discharge(u_obs, sp, lo, hi)
        rows.append((lab, m, o, m / o if o > 1e-9 else float("nan")))
    return {"path": path, "q_model": q_m, "q_obs": q_o,
            "ratio": q_m / q_o if q_o > 1e-9 else float("nan"), "bands": rows}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("maps", nargs="+")
    args = ap.parse_args()
    results = [score(p) for p in args.maps]
    print("\n" + "=" * 72)
    for r in results:
        print(f"\n{os.path.basename(r['path'])}")
        print(f"  discharge  model {r['q_model']:8.0f}   with u_obs "
              f"{r['q_obs']:8.0f} Gt/yr   ratio {r['ratio']:.2f}"
              f"   (observed {OBSERVED_DISCHARGE[0]:.0f} +/- "
              f"{OBSERVED_DISCHARGE[1]:.0f})")
        print("    discharge by observed speed of the source cell [Gt/yr]")
        print(f"    {'speed [m/yr]':>14s} {'model':>9s} {'u_obs':>9s} {'ratio':>7s}")
        for lab, m, o, ratio in r["bands"]:
            print(f"    {lab:>14s} {m:9.0f} {o:9.0f} {ratio:7.2f}")
    if len(results) > 1:
        print("\n  band ratios side by side")
        names = [os.path.basename(r["path"]) for r in results]
        print("    " + " ".join(f"{n[:18]:>18s}" for n in ["speed"] + names))
        for i, lab in enumerate(BAND_LABELS):
            cells = " ".join(f"{r['bands'][i][3]:18.2f}" for r in results)
            print(f"    {lab:>18s} {cells}")
        print("    " + f"{'overall':>18s} "
              + " ".join(f"{r['ratio']:18.2f}" for r in results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
