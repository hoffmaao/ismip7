#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 7: SSP5-8.5 with CESM2-WACCM (2015-2300).

Usage:
    mpiexec -n 12 python scripts/projections/ssp585_cesm_waccm.py
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
# dt=0.1 per the May 2026 dt sweep: captures shelf melt faithfully without
# oversampling the post-continuation transient.
DT = float(os.environ.get("ISMIP7_DT", "0.1"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

ESM = "CESM2-WACCM"
SSP = "ssp585"

_RESULTS = os.path.join(_PROJECT, "antarctica", "results")


def _find_k_npz(lc):
    r"""Calibrated per-basin K npz: prefer this mesh's calibration, fall back
    to the 2500 m one (per-basin K is 16 scalars remapped through the IMBIE2
    8 km grid, so it is mesh-independent; K* moved 3% between 2 and 2.5 km)."""
    override = os.environ.get("ISMIP7_K_PER_BASIN_NPZ")
    candidates = [override] if override else [
        os.path.join(_RESULTS, f"calibrated_K_per_basin_{lc}.npz"),
        os.path.join(_RESULTS, "calibrated_K_per_basin_2500.npz"),
    ]
    for c in candidates:
        if c and os.path.exists(c):
            return c
    return None


def main():
    restart = os.environ.get(
        "ISMIP7_RESTART",
        os.path.join(RESULTS_DIR, f"hist_cesm2_waccm_{lc}_final.h5"),
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

    K_npz = _find_k_npz(lc)
    K_melt = float(os.environ.get("ISMIP7_K_MELT", "1.15e-4"))
    if K_npz is not None:
        PETSc.Sys.Print(f"  Ocean melt: calibrated per-basin K from {K_npz}")
    else:
        PETSc.Sys.Print(
            f"  WARNING: no per-basin K calibration found; "
            f"scalar K={K_melt:.2e} everywhere"
        )
    # Full acabf (not the 1960-1989-referenced anomaly): the CTRL uses the
    # 2015-2029 acabf climatology, so projection minus CTRL is the forced
    # signal relative to the 2015-centered climate with no re-referencing.
    callback = make_forcing_callback(atm=atm, ocean=ocean, fracture=fracture,
                                     K=K_melt, K_per_basin_npz=K_npz,
                                     smb_anomaly=False)

    PETSc.Sys.Print(f"\nCore Experiment 7: {SSP} / {ESM}")
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
