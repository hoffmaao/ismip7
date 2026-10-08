#!/usr/bin/env python3
"""Where a Greene et al. (2022) ice mask and BedMachine disagree (issue #167).

The _front<year> meshes follow the mask's marine front, and the
greene<year> raster sampling takes front-cell thickness from BedMachine. This
counts, on the mask's 500 m pixels, the ice the mask holds where BedMachine
has no thickness (by BedMachine's mask there, and how much of it Bedmap3
covers) and the ice BedMachine holds outside the mask (by its mask, with its
volume), and describes the mask's marine front: its length and how much of
it BedMachine holds ice on.

    python antarctica/scripts/front_mask_census.py --year 2015 [--out census.json]

BedMachine comes from <ISMIP7_OBS_DATA_ROOT>/bedmachine, the mask from
obs_icemask.icemask_tif (cut from the MIPkit on first use), Bedmap3 from the
MIPkit when one is on this machine.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from icepack2_tools.obs_icemask import (  # noqa: E402
    ICE, LAND, MARINE, _aligned_window, icemask_source, icemask_tif, read_classes,
)
from icepack2_tools.runconfig import obs_data_root  # noqa: E402

RHO_I = 917.0


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--year", type=int, default=2015)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    import rasterio

    bm_fn = sorted(glob.glob(os.path.join(obs_data_root(), "bedmachine", "*.nc")))[0]
    tif = icemask_tif(args.year)
    kit, _ = icemask_source(args.year)
    classes, tr = read_classes(tif, bm_fn)
    bounds = rasterio.transform.array_bounds(*classes.shape, tr)
    window_bounds = (bounds[0], bounds[1], bounds[2], bounds[3])
    with rasterio.open(f"netcdf:{bm_fn}:thickness") as src:
        H = src.read(1, window=_aligned_window(src, window_bounds)).astype("f4")
    with rasterio.open(f"netcdf:{bm_fn}:mask") as src:
        bm_mask = src.read(1, window=_aligned_window(src, window_bounds))
    px_km2 = abs(tr.a * tr.e) / 1e6
    ice = classes == ICE
    held = H > 0

    def km2(sel):
        return round(float(sel.sum()) * px_km2, 1)

    def gt(sel):
        return round(float(H[sel].astype("f8").sum()) * px_km2 * 1e6 * RHO_I / 1e12, 1)

    out = {"mask": os.path.basename(tif), "bedmachine": os.path.basename(bm_fn),
           "pixel_km2": px_km2,
           "mask_ice_km2": km2(ice),
           "bedmachine_ice_km2": km2((bm_mask >= 2) & (bm_mask <= 4)),
           "mask_ice_bedmachine_empty_km2": {}, "bedmachine_ice_outside_mask": {}}
    empty = ice & ~held
    for name, v in (("ocean", 0), ("land", 1), ("grounded", 2), ("floating", 3), ("lake", 4)):
        out["mask_ice_bedmachine_empty_km2"][name] = km2(empty & (bm_mask == v))
    outside = ~ice & held
    for name, v in (("grounded", 2), ("floating", 3), ("lake", 4)):
        sel = outside & (bm_mask == v)
        out["bedmachine_ice_outside_mask"][name] = {"km2": km2(sel), "gt": gt(sel)}
    out["bedmachine_ice_outside_mask"]["marine_km2"] = km2(outside & (classes == MARINE))
    out["bedmachine_ice_outside_mask"]["land_km2"] = km2(outside & (classes == LAND))

    bm3 = None
    if kit is not None:
        import netCDF4
        with netCDF4.Dataset(kit) as d:
            if "thickness_bedmap3" in d.variables:
                v = d.variables["thickness_bedmap3"]
                v.set_auto_mask(False)
                bm3 = np.asarray(v[:], dtype="f4")
                if np.asarray(d.variables["y"][:2]).tolist() != sorted(
                        np.asarray(d.variables["y"][:2]).tolist(), reverse=True):
                    bm3 = bm3[::-1]
        if bm3 is not None and bm3.shape == classes.shape:
            covered = empty & np.isfinite(bm3) & (bm3 > 0)
            out["mask_ice_bedmachine_empty_bedmap3_km2"] = km2(covered)
            out["bedmap3_thickness_there_m"] = (
                round(float(np.median(bm3[covered])), 1) if covered.any() else None)

    front = np.zeros_like(ice)
    marine = classes == MARINE
    for di, dj in ((0, 1), (0, -1), (1, 0), (-1, 0)):
        nb = np.zeros_like(marine)
        src_r = slice(max(di, 0), marine.shape[0] + min(di, 0))
        dst_r = slice(max(-di, 0), marine.shape[0] + min(-di, 0))
        src_c = slice(max(dj, 0), marine.shape[1] + min(dj, 0))
        dst_c = slice(max(-dj, 0), marine.shape[1] + min(-dj, 0))
        nb[dst_r, dst_c] = marine[src_r, src_c]
        front |= ice & nb
    hf = H[front]
    out["marine_front"] = {
        "pixels": int(front.sum()),
        "bedmachine_holds_ice": round(float((hf > 0).mean()), 4),
        "thickness_held_mean_m": round(float(hf[hf > 0].mean()), 1),
        "thickness_held_median_m": round(float(np.median(hf[hf > 0])), 1),
        "by_bedmachine_mask": {name: int((front & (bm_mask == v)).sum())
                               for name, v in (("ocean", 0), ("land", 1), ("grounded", 2),
                                               ("floating", 3))},
    }
    # How far the two fronts sit apart: each pixel of mask ice over
    # BedMachine water from the nearest pixel BedMachine holds ice on, and
    # each pixel of BedMachine floating ice outside the mask from the
    # nearest mask ice, binned by area.
    from scipy.ndimage import distance_transform_edt
    with rasterio.open(f"netcdf:{bm_fn}:bed") as src:
        bed = src.read(1, window=_aligned_window(src, window_bounds))
    from icepack2_tools.obs_icemask import classify
    over_water = classify(np.zeros_like(ice), bm_mask, bed) == MARINE
    del bed
    px_km = abs(tr.a) / 1e3
    edges_km = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, np.inf]

    def binned(dist_km, sel):
        h, _ = np.histogram(dist_km[sel], bins=edges_km)
        return {f"{a:g}-{b:g}": round(float(c) * px_km2, 1)
                for a, b, c in zip(edges_km[:-1], edges_km[1:], h)}

    d_held = distance_transform_edt(~held) * px_km
    advanced = ice & ~held & over_water
    out["mask_ice_over_bedmachine_water_km2_by_km_to_bedmachine_ice"] = binned(d_held, advanced)
    out["marine_front_pixels_without_bedmachine_ice_km2_by_km_to_bedmachine_ice"] = binned(
        d_held, front & ~held)
    del d_held
    d_ice = distance_transform_edt(~ice) * px_km
    out["bedmachine_floating_outside_mask_km2_by_km_to_mask_ice"] = binned(
        d_ice, outside & (bm_mask == 3))
    del d_ice

    text = json.dumps(out, indent=1)
    print(text)
    if args.out:
        with open(args.out, "w") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    main()
