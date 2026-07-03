#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 7: SSP5-8.5 with CESM2-WACCM (2015-2300).

Usage:
    mpiexec -n 12 python scripts/projections/ssp585_cesm_waccm.py
"""

import os, sys
import numpy as np

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
    load_racmo_smb_climatology,
)

T_START = 2015.0
T_END = float(os.environ.get("ISMIP7_T_END", "2300"))
# dt=0.1 per the May 2026 dt sweep: captures shelf melt faithfully without
# oversampling the post-continuation transient.
DT = float(os.environ.get("ISMIP7_DT", "0.1"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

ESM = "CESM2-WACCM"
SSP = "ssp585"

# Reference-climate window — must mirror control/run.py so that
# projection minus CTRL is the forced signal.
CLIM_START = int(os.environ.get("ISMIP7_CLIM_START", "2000"))
CLIM_END = int(os.environ.get("ISMIP7_CLIM_END", "2029"))

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
    # SMB: RACMO climatology baseline + re-referenced ISMIP7 anomaly,
    # mirroring the CTRL. aSMB is wrt the ESM's 1960-1989 climatology, so
    # subtract its mean over the CTRL reference window and add RACMO:
    #   SMB(t) = RACMO_clim + aSMB(t) - mean(aSMB | window)
    # Fallback without RACMO: full acabf(t) — the same construction with
    # the ESM's own acabf climatology as the baseline.
    smb_anomaly, smb_baseline = False, None
    if atm is not None:
        mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
        mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]
        ref_years = [y for y in atm.available_years("acabf-anomaly")
                     if CLIM_START <= y <= CLIM_END]
        try:
            if not ref_years:
                raise FileNotFoundError(
                    f"no acabf-anomaly years in {CLIM_START}-{CLIM_END}"
                )
            racmo = load_racmo_smb_climatology(ctx["Q"], CLIM_START, CLIM_END)
            ref = np.zeros(len(mesh_x))
            for yr in ref_years:
                ref += atm.get_smb(yr, mesh_x, mesh_y, anomaly=True)
            ref /= len(ref_years)
            smb_anomaly = True
            smb_baseline = racmo.dat.data_ro - ref
            PETSc.Sys.Print(
                f"  SMB: RACMO2.4p1 baseline + aSMB re-referenced to "
                f"{ref_years[0]}-{ref_years[-1]} ({len(ref_years)} yr)"
            )
        except FileNotFoundError as e:
            PETSc.Sys.Print(f"  SMB: no RACMO baseline ({e}); "
                            f"forcing with full acabf(t)")
    callback = make_forcing_callback(atm=atm, ocean=ocean, fracture=fracture,
                                     K=K_melt, K_per_basin_npz=K_npz,
                                     smb_anomaly=smb_anomaly,
                                     smb_baseline=smb_baseline)

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
