#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 11: OCX observationally constrained (1990-2025).

Usage:
    mpiexec -n 12 python scripts/projections/ocx.py
"""

import os, sys

_PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)
))))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from simulation import (
    setup_model, run_simulation, RESULTS_DIR, PETSc, lc,
)
from icepack2_tools.forcing import (
    ISMIP7Atmosphere, ISMIP7Ocean, make_forcing_callback,
)

T_START = float(os.environ.get("ISMIP7_T_START", "1990"))
T_END = float(os.environ.get("ISMIP7_T_END", "2025"))
DT = float(os.environ.get("ISMIP7_DT", "1.0"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "5"))


def main():
    ctx = setup_model()

    atm = ISMIP7Atmosphere(scenario="ocx")
    years = atm.available_years()
    if years:
        PETSc.Sys.Print(f"  Atmosphere forcing: {len(years)} years ({years[0]}-{years[-1]})")
    else:
        PETSc.Sys.Print(f"  WARNING: no atmosphere data found, using zero SMB anomaly")
        atm = None

    ocean = ISMIP7Ocean(scenario="ocx")
    callback = make_forcing_callback(atm=atm, ocean=ocean)

    PETSc.Sys.Print(f"\nCore Experiment 11: OCX (observationally constrained)")
    PETSc.Sys.Print(f"  Period: {T_START}-{T_END}")

    run_simulation(
        ctx,
        experiment_name="ocx",
        t_start=T_START,
        t_end=T_END,
        dt=DT,
        output_interval=OUTPUT_INTERVAL,
        forcing_callback=callback,
    )

    if ocean is not None:
        ocean.close()


if __name__ == "__main__":
    main()
