#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 6: SSP1-2.6 with MRI-ESM2-0 (2015-2300).

Usage:
    mpiexec -n 12 python scripts/projections/ssp126_mri_esm2.py
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
    ISMIP7Atmosphere, ISMIP7Ocean, ISMIP7Fracture, make_forcing_callback,
)

T_START = 2015.0
T_END = float(os.environ.get("ISMIP7_T_END", "2300"))
DT = float(os.environ.get("ISMIP7_DT", "1.0"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

ESM = "MRI-ESM2-0"
SSP = "ssp126"


def main():
    restart = os.environ.get(
        "ISMIP7_RESTART",
        os.path.join(RESULTS_DIR, f"hist_mri_esm2_0_{lc}_final.h5"),
    )
    if os.path.exists(restart):
        ctx = setup_model(restart_from=restart)
    else:
        PETSc.Sys.Print(f"  No historical restart found at {restart}")
        PETSc.Sys.Print(f"  Starting from BedMachine initial state")
        ctx = setup_model()

    atm = ISMIP7Atmosphere(esm=ESM, scenario=SSP)
    years = atm.available_years()
    if years:
        PETSc.Sys.Print(f"  Atmosphere forcing: {len(years)} years ({years[0]}-{years[-1]})")
    else:
        PETSc.Sys.Print(f"  WARNING: no atmosphere data found, using zero SMB anomaly")
        atm = None

    ocean = ISMIP7Ocean(esm=ESM, scenario=SSP)
    fracture = ISMIP7Fracture(esm=ESM, scenario=SSP)
    fracture.load()

    callback = make_forcing_callback(atm=atm, ocean=ocean, fracture=fracture)

    PETSc.Sys.Print(f"\nCore Experiment 6: {SSP} / {ESM}")
    PETSc.Sys.Print(f"  Period: {T_START}-{T_END}")

    run_simulation(
        ctx,
        experiment_name=f"{SSP}_{ESM.lower().replace('-', '_')}",
        t_start=T_START,
        t_end=T_END,
        dt=DT,
        output_interval=OUTPUT_INTERVAL,
        forcing_callback=callback,
    )

    if ocean is not None:
        ocean.close()
    if fracture is not None:
        fracture.close()


if __name__ == "__main__":
    main()
