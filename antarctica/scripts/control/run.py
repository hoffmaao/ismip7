#!/usr/bin/env python3
r"""ISMIP7 Core Experiments 9/10: CTRL2015 -- constant 2015 climate.

Atmosphere: 2000-2029 SMB climatology held fixed.
Ocean:      OI climatology TF + so held fixed; melt is recomputed each
            step from the evolving geometry using the Burgard quadratic-
            mixed-slope formula with per-basin calibrated K (see
            `antarctica/scripts/calibrate_melt.py`).

Usage:
    mpiexec -n 12 python scripts/control/run.py
    ISMIP7_ESM=MRI-ESM2-0 mpiexec -n 12 python scripts/control/run.py
"""

import os, sys, glob
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)
))))
sys.path.insert(0, _PROJECT)

from simulation import setup_model, run_simulation, PETSc, lc
from icepack2_tools.forcing import (
    ISMIP7Atmosphere,
    compute_sin_alpha,
    quadratic_mixed_slope,
    load_K_per_basin,
    _RHO_ICE, _RHO_WATER,
)

T_START = 2015.0
T_END = float(os.environ.get("ISMIP7_T_END", "2300"))
DT = float(os.environ.get("ISMIP7_DT", "1.0"))
OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))

ESM = os.environ.get("ISMIP7_ESM", "CESM2-WACCM")
# Reference-climate window for the constant SMB. Default 2015-2029: the
# local data tree only has the projection scenario (2015+), and the ssp585
# acabf-anomaly is referenced to 1960-1989, so an anomaly-based baseline is
# not constructible. Widen (e.g. 2000-2029) once historical acabf is
# downloaded — compute_climatology pools historical + projection years.
CLIM_START = int(os.environ.get("ISMIP7_CLIM_START", "2015"))
CLIM_END = int(os.environ.get("ISMIP7_CLIM_END", "2029"))
CLIM_SCENARIO = os.environ.get("ISMIP7_CLIM_SCENARIO", "ssp585")

DATA_ROOT = os.environ.get(
    "ISMIP7_DATA_ROOT", os.path.join(_PROJECT, "ISMIP7", "AIS")
)
CLIM_TF = os.path.join(DATA_ROOT, "meltMIP", "OI_Climatology_ismip8km_60m_tf_extrap.nc")
CLIM_SO = os.path.join(DATA_ROOT, "meltMIP", "OI_Climatology_ismip8km_60m_so_extrap.nc")
# Per-basin K: prefer this mesh's calibration, else the 2500 m one (16
# basin scalars remapped through the IMBIE2 8 km grid — mesh-independent).
_K_LC_NPZ = os.path.join(_PROJECT, "antarctica", "results",
                         f"calibrated_K_per_basin_{lc}.npz")
_K_2500_NPZ = os.path.join(_PROJECT, "antarctica", "results",
                           "calibrated_K_per_basin_2500.npz")
K_NPZ = os.environ.get(
    "ISMIP7_K_PER_BASIN_NPZ",
    _K_LC_NPZ if os.path.exists(_K_LC_NPZ) else _K_2500_NPZ,
)


def compute_climatology(atms, mesh_x, mesh_y):
    r"""Mean full-field (acabf) SMB over [CLIM_START, CLIM_END], pooling
    years across the given atmospheres (historical + projection scenario).

    Full field only: a mean of `acabf-anomaly` is an anomaly wrt the ESM's
    1960-1989 climatology, not a baseline, so no anomaly fallback exists.
    """
    smb_sum = np.zeros(len(mesh_x))
    n = 0
    for atm in atms:
        years = [y for y in atm.available_years("acabf")
                 if CLIM_START <= y <= CLIM_END]
        for yr in years:
            smb_sum += atm.get_smb(yr, mesh_x, mesh_y, anomaly=False)
        n += len(years)
        if years:
            PETSc.Sys.Print(
                f"  {atm.scenario}: {len(years)} acabf years "
                f"({years[0]}-{years[-1]}) in climatology window"
            )
    if n == 0:
        return None
    return smb_sum / n


def _build_climatology_interpolators():
    r"""Load OI climatology TF and so into RegularGridInterpolator (z, y, x)."""
    interps = {}
    for path, var in [(CLIM_TF, "tf"), (CLIM_SO, "so")]:
        ds = xr.open_dataset(path)
        da = ds[var]
        zdim = [d for d in da.dims if d.lower() in ("z", "depth", "lev")][0]
        za = ds[zdim].values
        ya = ds["y"].values
        xa = ds["x"].values
        data = da.transpose(zdim, "y", "x").values.astype(np.float32)
        if za[0] > za[-1]:
            za = za[::-1]; data = data[::-1, :, :]
        if ya[0] > ya[-1]:
            ya = ya[::-1]; data = data[:, ::-1, :]
        if xa[0] > xa[-1]:
            xa = xa[::-1]; data = data[:, :, ::-1]
        data = np.nan_to_num(data, nan=0.0)
        interps[var] = (
            RegularGridInterpolator(
                (za, ya, xa), data,
                method="nearest", bounds_error=False, fill_value=0.0,
            ),
            za,
        )
        ds.close()
    return interps


def make_ctrl_ocean_callback(K_field):
    r"""Build a CTRL2015 ocean-melt callback: constant climatology TF/so,
    evolving geometry, per-node K from calibration."""
    PETSc.Sys.Print("  Building OI-climatology interpolators...")
    interps = _build_climatology_interpolators()

    def callback(ctx, t_yr):
        mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
        mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]
        h = ctx["h"].dat.data_ro
        b = ctx["b"].dat.data_ro
        s = ctx["s"].dat.data_ro
        draft = np.minimum(s - h, 0.0)

        tf_interp, za_tf = interps["tf"]
        so_interp, za_so = interps["so"]
        d_tf = np.clip(draft, za_tf[0], za_tf[-1])
        d_so = np.clip(draft, za_so[0], za_so[-1])
        tf = tf_interp(np.column_stack([d_tf, mesh_y, mesh_x]))
        sal = so_interp(np.column_stack([d_so, mesh_y, mesh_x]))
        sin_a = compute_sin_alpha(ctx)

        melt = quadratic_mixed_slope(tf, sal, sin_a, K=K_field)

        haf = s - (b + (_RHO_WATER / _RHO_ICE) * np.maximum(-b, 0.0))
        floating = haf <= 0
        ctx["ocean_melt"].dat.data[:] = np.where(floating, melt, 0.0)

    return callback


def main():
    ctx = setup_model()

    mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
    mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]

    # Atmosphere: fixed SMB climatology, pooled over historical + the
    # projection scenario (whichever years exist locally in the window).
    atms = [
        ISMIP7Atmosphere(esm=ESM, scenario="historical"),
        ISMIP7Atmosphere(esm=ESM, scenario=CLIM_SCENARIO),
    ]
    clim_smb = compute_climatology(atms, mesh_x, mesh_y)
    if clim_smb is None:
        if os.environ.get("ISMIP7_ALLOW_ZERO_SMB"):
            PETSc.Sys.Print("  WARNING: no climatology data, using zero SMB")
            ctx["accum"].assign(0.0)
        else:
            raise FileNotFoundError(
                f"No acabf data for {ESM} in {CLIM_START}-{CLIM_END} "
                f"(historical or {CLIM_SCENARIO}). A zero-SMB control is "
                f"almost certainly not what you want; set "
                f"ISMIP7_ALLOW_ZERO_SMB=1 to force it."
            )
    else:
        ctx["accum"].dat.data[:] = clim_smb
        PETSc.Sys.Print(f"  Climatological SMB: mean={clim_smb.mean():.4f} m/yr")

    # Ocean: fixed climatology TF/so + per-basin calibrated K
    if not os.path.exists(K_NPZ):
        raise FileNotFoundError(
            f"Per-basin K calibration not found at {K_NPZ}. "
            f"Run antarctica/scripts/calibrate_melt.py first."
        )
    PETSc.Sys.Print(f"  Loading per-basin K from: {K_NPZ}")
    K_field = load_K_per_basin(K_NPZ, mesh_x, mesh_y, fill=0.0)
    PETSc.Sys.Print(
        f"  K field: nonzero={int((K_field>0).sum())}/{len(K_field)}  "
        f"med={np.median(K_field[K_field>0]) if (K_field>0).any() else 0:.2e}"
    )

    callback = make_ctrl_ocean_callback(K_field)

    PETSc.Sys.Print(f"\nControl experiment: {ESM}")
    PETSc.Sys.Print(f"  Period: {T_START}-{T_END}")
    PETSc.Sys.Print(f"  Constant {CLIM_START}-{CLIM_END} SMB climatology")
    PETSc.Sys.Print(f"  Constant OI ocean climatology + per-basin K")

    esm_tag = ESM.lower().replace("-", "_")
    run_simulation(
        ctx,
        experiment_name=f"ctrl2015_{esm_tag}",
        t_start=T_START,
        t_end=T_END,
        dt=DT,
        output_interval=OUTPUT_INTERVAL,
        forcing_callback=callback,
    )


if __name__ == "__main__":
    main()
