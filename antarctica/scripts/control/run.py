#!/usr/bin/env python3
r"""ISMIP7 Core Experiments 9/10: CTRL2015 -- constant 2015 climate.

Usage:
    mpiexec -n 12 python scripts/control/run.py
    ISMIP7_ESM=MRI-ESM2-0 mpiexec -n 12 python scripts/control/run.py
"""

import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)
))))
sys.path.insert(0, _PROJECT)

from firedrake import assemble, dx, Constant
from simulation import setup_model, run_simulation, PETSc
from icepack2_tools.forcing import (
    ISMIP7Atmosphere, smb_kgm2s_to_myr, load_racmo_smb_climatology,
)

T_START = 2015.0
T_END = float(os.environ.get("ISMIP7_T_END", "2300"))
DT = float(os.environ.get("ISMIP7_DT", "1.0"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

ESM = os.environ.get("ISMIP7_ESM", "CESM2-WACCM")
CLIM_START = 2000
CLIM_END = 2029


def area_weighted_mean(field, mesh):
    r"""Domain-area-weighted mean of a field (mesh refinement varies spatially)."""
    return assemble(field * dx) / assemble(Constant(1.0) * dx(domain=mesh))


def compute_climatology(atm, mesh_x, mesh_y):
    r"""Compute 2000-2029 mean SMB from available acabf files."""
    years = atm.available_years("acabf")
    clim_years = [y for y in years if CLIM_START <= y <= CLIM_END]

    if not clim_years:
        PETSc.Sys.Print(f"  No acabf data for {CLIM_START}-{CLIM_END}")
        years_anom = atm.available_years("acabf-anomaly")
        clim_years_anom = [y for y in years_anom if CLIM_START <= y <= CLIM_END]
        if clim_years_anom:
            PETSc.Sys.Print(f"  Using acabf-anomaly mean over {clim_years_anom[0]}-{clim_years_anom[-1]}")
            smb_sum = np.zeros(len(mesh_x))
            for yr in clim_years_anom:
                smb_sum += atm.get_smb(yr, mesh_x, mesh_y, anomaly=True)
            return smb_sum / len(clim_years_anom)
        return None

    PETSc.Sys.Print(f"  Computing {CLIM_START}-{CLIM_END} climatology from {len(clim_years)} years")
    smb_sum = np.zeros(len(mesh_x))
    for yr in clim_years:
        smb_sum += atm.get_smb(yr, mesh_x, mesh_y, anomaly=False)
    return smb_sum / len(clim_years)


def main():
    ctx = setup_model()

    atm = ISMIP7Atmosphere(esm=ESM, scenario="historical")
    mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
    mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]

    try:
        ctx["accum"].assign(
            load_racmo_smb_climatology(ctx["Q"], CLIM_START, CLIM_END)
        )
        mean_smb = area_weighted_mean(ctx["accum"], ctx["mesh"])
        PETSc.Sys.Print(
            f"  Climatological SMB from RACMO2.4p1 ({CLIM_START}-{CLIM_END}): "
            f"area-weighted mean={mean_smb:.4f} m/yr"
        )
    except FileNotFoundError:
        PETSc.Sys.Print("  No RACMO data; falling back to ISMIP7 acabf climatology")
        clim_smb = compute_climatology(atm, mesh_x, mesh_y)
        if clim_smb is not None:
            ctx["accum"].dat.data[:] = clim_smb
            mean_smb = area_weighted_mean(ctx["accum"], ctx["mesh"])
            PETSc.Sys.Print(f"  Climatological SMB: area-weighted mean={mean_smb:.4f} m/yr")
        else:
            PETSc.Sys.Print("  WARNING: no climatology data, using zero SMB")
            ctx["accum"].assign(0.0)

    PETSc.Sys.Print(f"\nControl experiment: {ESM}")
    PETSc.Sys.Print(f"  Period: {T_START}-{T_END}")
    PETSc.Sys.Print(f"  Constant {CLIM_START}-{CLIM_END} climatology")

    esm_tag = ESM.lower().replace("-", "_")
    run_simulation(
        ctx,
        experiment_name=f"ctrl2015_{esm_tag}",
        t_start=T_START,
        t_end=T_END,
        dt=DT,
        output_interval=OUTPUT_INTERVAL,
        forcing_callback=None,
    )


if __name__ == "__main__":
    main()
