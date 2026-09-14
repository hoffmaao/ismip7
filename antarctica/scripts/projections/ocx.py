#!/usr/bin/env python3
r"""ISMIP7 Core Experiment 11: OCX observationally constrained (1979-2025).

Fully observation-forced run over the satellite era, for validating the
initialized model against the observed record:

    SMB:   RACMO2.4p1 actual-year fields (1979-2023; end years held at the
           last available RACMO year), NOT an ESM.
    Ocean: constant OI-climatology TF/so at draft with the calibrated
           per-basin K (the same forcing the CTRL uses).

If an `ocx` scenario tree exists under ISMIP7/AIS (protocol-provided
time-varying obs forcing), the atmosphere/ocean readers pick it up
instead. The period is the OCX forcing span (the OCX atmosphere tree runs
1979-2025). Note the initial state is the ~2015 BedMachine/MAP geometry, so
a 1979 start is anachronistic by construction: treat the early years as
relaxation and the 2000s-2025 as the validation window.

Usage:
    mpiexec -n 24 python scripts/projections/ocx.py
"""
import os, sys

_SCRIPTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from simulation import (setup_model, run_simulation, latest_checkpoint,
                        auto_resume, PETSc)
from experiment import find_k_npz
from icepack2_tools.forcing import (
    ISMIP7Atmosphere, ISMIP7Ocean, make_forcing_callback,
    make_climatology_ocean_callback, load_racmo_smb_climatology,
    load_K_per_basin, forcing_coords, reject_collapse_mask, forcing_year,
)

T_START = float(os.environ.get("ISMIP7_T_START", "1979"))
# 1 January of the year AFTER the last one covered, the convention every core
# driver uses: years 1979 through 2025 run and 2025 is the last banked year.
T_END = float(os.environ.get("ISMIP7_T_END", "2026"))
DT = float(os.environ.get("ISMIP7_DT", "0.1"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))
RACMO_LAST = 2023  # smbgl_monthlyS_ANT11_RACMO2.4p1_ERA5_197901_202312


def main():
    experiment_name = "ocx" + (f"_{os.environ.get('ISMIP7_RUN_TAG', '')}"
                               if os.environ.get("ISMIP7_RUN_TAG") else "")
    # Explicit restart, else unattended auto-resume from this experiment's own
    # newest checkpoint (ISMIP7_AUTO_RESUME=1), the same lookup the control
    # driver does. A chained batch job depends on it: without it every link
    # cold-starts and the chain never advances.
    restart = os.environ.get("ISMIP7_RESTART")
    if restart is None and auto_resume():
        restart = latest_checkpoint(experiment_name)
        PETSc.Sys.Print(
            f"Auto-resume: {restart}" if restart
            else "Auto-resume: no prior checkpoint"
        )
    ctx = setup_model(restart_from=restart)
    # Sample forcing at the geometry dofs, not the mesh vertices: under
    # DG0 geometry those are cell centroids (see forcing.forcing_coords).
    mesh_x, mesh_y = forcing_coords(ctx)

    # Per-basin K (with 2500 m fallback) + optional global scale.
    K_npz = find_k_npz()
    if K_npz is None:
        raise FileNotFoundError(
            "OCX needs the calibrated per-basin K "
            "(antarctica/scripts/calibrate_melt.py)."
        )
    K_field = load_K_per_basin(K_npz, mesh_x, mesh_y, fill=0.0)
    K_scale = float(os.environ.get("ISMIP7_K_SCALE", "1.0"))
    if K_scale != 1.0:
        K_field = K_field * K_scale
        PETSc.Sys.Print(f"  K scaled by ISMIP7_K_SCALE={K_scale:.3f}")
    PETSc.Sys.Print(f"  Ocean melt: OI climatology + per-basin K ({K_npz})")

    reject_collapse_mask("the OCX experiment")

    atm = ISMIP7Atmosphere(scenario="ocx")
    if atm.available_years():
        PETSc.Sys.Print("  Atmosphere: protocol ocx tree")
        ocean = ISMIP7Ocean(scenario="ocx")
        callback = make_forcing_callback(
            atm=atm, ocean=ocean, K_per_basin_npz=K_npz, smb_anomaly=False,
        )
    else:
        PETSc.Sys.Print("  Atmosphere: RACMO2.4p1 actual-year SMB (no ocx tree)")
        ocean = None
        oi_melt = make_climatology_ocean_callback(K_field)
        racmo_cache = {}

        def callback(ctx_, t_yr):
            yr = min(forcing_year(t_yr), RACMO_LAST)
            if yr not in racmo_cache:
                racmo_cache[yr] = load_racmo_smb_climatology(
                    ctx_["Q_g"], yr, yr
                ).dat.data_ro.copy()
                while len(racmo_cache) > 3:
                    racmo_cache.pop(next(iter(racmo_cache)))
            ctx_["accum"].dat.data[:] = racmo_cache[yr]
            oi_melt(ctx_, t_yr)

    PETSc.Sys.Print("\nCore Experiment 11: OCX (observationally constrained)")
    PETSc.Sys.Print(f"  Period: {T_START}-{T_END}, dt={DT}")

    run_simulation(
        ctx,
        experiment_name=experiment_name,
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
