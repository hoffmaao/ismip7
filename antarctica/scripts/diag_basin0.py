#!/usr/bin/env python3
r"""Quick diagnostic: where is IMBIE2 basin 0?

Plots, on the ISMIP7 8 km grid:
 - basin 0 cells (vs the named basins 1-15)
 - BedMachine v4.1 floating-ice mask remapped to 8 km
 - the intersection (basin 0 AND floating) - these are the cells that
   pollute the calibration when they should belong to a named basin

Saves antarctica/results/diag_basin0.png.
"""

import os, sys
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(os.path.dirname(_HERE))

DATA = os.path.join(_PROJECT, "ISMIP7", "AIS")
IMBIE2 = os.path.join(DATA, "parameterisations", "ocean", "imbie2",
                      "basin_numbers_ismip8km_v2.nc")
import glob
BM = glob.glob(os.path.join(_PROJECT, "antarctica", "data",
                            "bedmachine", "*.nc"))[0]
OUT = os.path.join(_PROJECT, "antarctica", "results", "diag_basin0.png")
os.makedirs(os.path.dirname(OUT), exist_ok=True)


def main():
    # IMBIE2 basins on 8km grid
    bds = xr.open_dataset(IMBIE2)
    basin = bds["basinNumber"].values
    x8 = bds["x"].values
    y8 = bds["y"].values

    # BedMachine: load mask, remap to 8km via nearest-neighbour
    ds = xr.open_dataset(BM)
    bm_mask = ds["mask"].values  # 0=ocean, 1=land, 2=grounded, 3=floating, 4=Vostok
    bm_x = ds["x"].values
    bm_y = ds["y"].values
    if bm_y[0] > bm_y[-1]:
        bm_y = bm_y[::-1]; bm_mask = bm_mask[::-1, :]

    interp = RegularGridInterpolator(
        (bm_y, bm_x), bm_mask.astype(np.float32),
        method="nearest", bounds_error=False, fill_value=0.0,
    )
    Y8, X8 = np.meshgrid(y8, x8, indexing="ij")
    bm_on_8km = interp(np.column_stack([Y8.ravel(), X8.ravel()])).reshape(Y8.shape)
    floating_8km = (np.round(bm_on_8km).astype(int) == 3)

    # Cells: basin 0, basin 0 AND floating, named basins
    is_b0 = (basin == 0)
    is_b0_and_float = is_b0 & floating_8km
    is_named = (basin > 0)
    is_named_float = is_named & floating_8km

    n_b0_float = int(is_b0_and_float.sum())
    n_named_float = int(is_named_float.sum())
    print(f"  8km cells with basin=0: {int(is_b0.sum())}")
    print(f"  8km floating cells: {int(floating_8km.sum())}")
    print(f"  basin=0 AND floating: {n_b0_float}  "
          f"(area ~{n_b0_float*64} km^2)")
    print(f"  basin>0 AND floating: {n_named_float}  "
          f"(area ~{n_named_float*64} km^2)")
    print(f"  ratio basin0_float / total_float: "
          f"{n_b0_float / max(n_b0_float + n_named_float, 1):.4f}")

    # Per named basin counts
    print("  named-basin floating cell counts:")
    for b in range(1, 16):
        c = int(((basin == b) & floating_8km).sum())
        print(f"    basin {b:2d}: {c:6d} cells = ~{c*64:7d} km^2")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(16, 6),
                             constrained_layout=True)
    extent = [x8[0]/1e3, x8[-1]/1e3, y8[0]/1e3, y8[-1]/1e3]

    ax = axes[0]
    im = ax.imshow(basin, extent=extent, origin="lower",
                   cmap="tab20", vmin=0, vmax=20)
    ax.set_title("IMBIE2 basinNumber")
    plt.colorbar(im, ax=ax, fraction=0.046)
    ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]")

    ax = axes[1]
    cmap = matplotlib.colors.ListedColormap(
        ["lightgray", "lightblue", "lightgreen", "tomato", "purple"]
    )
    ax.imshow(np.round(bm_on_8km).astype(int), extent=extent, origin="lower",
              cmap=cmap, vmin=0, vmax=4)
    ax.set_title("BedMachine mask on 8km\n0=ocean 1=land 2=grounded 3=floating 4=lake")
    ax.set_xlabel("x [km]")

    ax = axes[2]
    # Show: red where basin==0 & floating; green where basin>0 & floating
    rgb = np.zeros((*basin.shape, 3), dtype=np.float32)
    rgb[..., 0] = is_b0_and_float.astype(np.float32)     # red
    rgb[..., 1] = is_named_float.astype(np.float32)      # green
    ax.imshow(rgb, extent=extent, origin="lower")
    ax.set_title(f"Floating cells: red=basin 0 ({n_b0_float})  "
                 f"green=basin 1-15 ({n_named_float})")
    ax.set_xlabel("x [km]")

    fig.suptitle("Basin-0 diagnostic — ISMIP7 8 km grid", fontsize=14)
    fig.savefig(OUT, dpi=120, bbox_inches="tight")
    print(f"  Saved: {OUT}")


if __name__ == "__main__":
    main()
