#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 2: Historical with MRI-ESM2-0 (>=1850 to 2014).

Usage:
    mpiexec -n 12 python scripts/historical/mri_esm2.py
"""

import os, sys

_PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)
))))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from simulation import setup_model, run_simulation, PETSc
from icepack2_tools.forcing import (
    ISMIP7Atmosphere, ISMIP7Ocean, make_forcing_callback,
)

T_START = float(os.environ.get("ISMIP7_T_START", "1850"))
T_END = 2014.0
DT = float(os.environ.get("ISMIP7_DT", "1.0"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

ESM = "MRI-ESM2-0"


def main():
    ctx = setup_model()

    atm = ISMIP7Atmosphere(esm=ESM, scenario="historical")
    years = atm.available_years()
    if years:
        PETSc.Sys.Print(f"  Atmosphere forcing: {len(years)} years ({years[0]}-{years[-1]})")
    else:
        PETSc.Sys.Print(f"  WARNING: no atmosphere data found, using zero SMB anomaly")
        atm = None

    ocean = ISMIP7Ocean(esm=ESM, scenario="historical")
    callback = make_forcing_callback(atm=atm, ocean=ocean)

    PETSc.Sys.Print(f"\nCore Experiment 2: Historical / {ESM}")
    PETSc.Sys.Print(f"  Period: {T_START}-{T_END}")

    run_simulation(
        ctx,
        experiment_name=f"hist_{ESM.lower().replace('-', '_')}",
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
