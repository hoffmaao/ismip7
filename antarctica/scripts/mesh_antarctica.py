#!/usr/bin/env python3
"""
Generate isotropic Antarctic mesh with Ua-style grounding zone refinement.

Refinement based on distance from the grounding line (Ua GLrange approach)
and strain-rate based sizing for outlet glaciers.

Usage:
    ISMIP7_BUFFER_M=20000 python scripts/mesh_antarctica.py
"""

import os, sys, glob
import numpy as np
import xarray as xr
import gmsh

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PROJECT = os.path.dirname(_ROOT)
sys.path.insert(0, _PROJECT)

from icepack2_tools.mesh import (
    load_bedmachine_mask,
    extract_ice_outline,
    classify_boundaries,
    load_velocity_for_sizing,
    build_gmsh_geometry,
    SUBSAMPLE,
)

DATA_DIR = os.path.join(_ROOT, "data")
MESH_DIR = os.path.join(_ROOT, "mesh")

# Ua-style distance bands: (distance_m, element_size_m)
GL_BANDS = [
    (5e3, 2500),
    (15e3, 5000),
    (30e3, 8000),
]
SHELF_SIZE = 5000
BUFFER_SIZE = 10000
COARSE = 64000


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def grounding_line_distance():
    """Distance from the grounding line (HAF = 0 boundary)."""
    from scipy.ndimage import distance_transform_edt, binary_dilation

    print("Computing GL distance...")
    fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    ds = xr.open_dataset(fn)
    thick = ds["thickness"].values[::SUBSAMPLE, ::SUBSAMPLE]
    bed = ds["bed"].values[::SUBSAMPLE, ::SUBSAMPLE]
    mask_bm = ds["mask"].values[::SUBSAMPLE, ::SUBSAMPLE]
    x = ds["x"].values[::SUBSAMPLE]
    y = ds["y"].values[::SUBSAMPLE]
    ds.close()
    dx = abs(float(x[1] - x[0]))

    rho_I, rho_W = 917.0, 1024.0
    s = np.maximum(bed + thick, (1 - rho_I / rho_W) * thick)
    s_float = bed + (rho_W / rho_I) * np.maximum(-bed, 0)
    haf = s - s_float

    ice = (mask_bm >= 2) & (mask_bm <= 4)
    grounded = ice & (haf > 0)
    floating = ice & (haf <= 0)

    struct = np.ones((3, 3), dtype=bool)
    gl_region = binary_dilation(grounded, struct, 2) & binary_dilation(floating, struct, 2) & ice

    dist = distance_transform_edt(~gl_region) * dx
    has_ice = (thick > 10) & ice
    is_float = mask_bm == 3

    coords = {"x": x.astype(float), "y": y.astype(float)}
    dist_da = xr.DataArray(dist, dims=["y", "x"], coords=coords)
    ice_da = xr.DataArray(has_ice.astype(float), dims=["y", "x"], coords=coords)
    float_da = xr.DataArray(is_float.astype(float), dims=["y", "x"], coords=coords)

    print(f"  GL region: {gl_region.sum()} pixels")
    return dist_da, ice_da, float_da


def calving_front_distance(boundaries, names):
    """Distance from calving boundary segments + BedMachine ice edge."""
    from scipy.ndimage import distance_transform_edt

    print("Computing calving front distance...")
    fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    ds = xr.open_dataset(fn)
    thick = ds["thickness"].values[::SUBSAMPLE, ::SUBSAMPLE]
    mask_bm = ds["mask"].values[::SUBSAMPLE, ::SUBSAMPLE]
    x = ds["x"].values[::SUBSAMPLE]
    y = ds["y"].values[::SUBSAMPLE]
    ds.close()
    dx = abs(float(x[1] - x[0]))

    has_ice = (thick > 10) & (mask_bm >= 2) & (mask_bm <= 4)
    dist_into = distance_transform_edt(has_ice) * dx
    ice_edge = has_ice & (dist_into < dx * 3)

    for seg_coords, name in zip(boundaries, names):
        if not name.startswith("Calving"):
            continue
        for pt in seg_coords:
            ix = np.argmin(np.abs(x - pt[0]))
            iy = np.argmin(np.abs(y - pt[1]))
            for di in range(-3, 4):
                for dj in range(-3, 4):
                    ii, jj = iy + di, ix + dj
                    if 0 <= ii < len(y) and 0 <= jj < len(x):
                        ice_edge[ii, jj] = True

    dist = distance_transform_edt(~ice_edge) * dx
    return xr.DataArray(dist, dims=["y", "x"],
                        coords={"x": x.astype(float), "y": y.astype(float)})


def main():
    os.makedirs(MESH_DIR, exist_ok=True)

    buffer_m = float(os.environ.get("ISMIP7_BUFFER_M", "20000"))
    if buffer_m == 0:
        os.environ["ISMIP7_BUFFER_M"] = "20000"
    print(f"Buffer: {buffer_m/1e3:.0f} km")

    mask, x, y = load_bedmachine_mask()
    outline = extract_ice_outline(mask, x, y)
    boundaries, names = classify_boundaries(outline, mask, x, y)
    refinement = load_velocity_for_sizing()

    gl_dist, ice_field, float_field = grounding_line_distance()
    cf_dist = calving_front_distance(boundaries, names)

    fn_base = os.path.join(MESH_DIR, f"antarctica_{COARSE}_{GL_BANDS[0][1]}")

    # Pass 1: raw mesh
    print(f"\nPass 1: raw mesh...")
    gmsh.initialize(sys.argv)
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add(fn_base + "_raw")

    build_gmsh_geometry(boundaries, names, 4000, COARSE)
    gmsh.model.mesh.generate(2)
    gmsh.write(fn_base + "_raw.msh")

    vtags, vxyz, _ = gmsh.model.mesh.getNodes()
    vxyz = vxyz.reshape((-1, 3))
    vmap = {int(j): i for i, j in enumerate(vtags)}
    tri_tags, tri_vtags = gmsh.model.mesh.getElementsByType(2)
    tri_vids = np.array([vmap[int(j)] for j in tri_vtags])
    triangles = tri_vids.reshape((tri_tags.shape[-1], -1))
    tri_centers = vxyz[triangles].mean(axis=1)

    print(f"  {len(vtags)} nodes, {len(tri_tags)} triangles")

    cx = xr.DataArray(tri_centers[:, 0], dims="tri")
    cy = xr.DataArray(tri_centers[:, 1], dims="tri")

    def interp(da):
        fx = da.coords[da.dims[-1]]
        fy = da.coords[da.dims[-2]]
        return np.nan_to_num(
            da.interp({fx.name: cx, fy.name: cy}, method="nearest").values.flatten(),
            nan=0.0,
        )

    ref_vals = np.maximum(interp(refinement), 1e-8)
    gl_d = interp(gl_dist)
    cf_d = interp(cf_dist)
    on_ice = interp(ice_field) > 0.5
    on_shelf = interp(float_field) > 0.5

    # Build size field
    # 1. Strain rate
    sr_floor = 8000
    size_sr = np.clip(sr_floor / ref_vals, sr_floor, COARSE)

    # 2. GL distance bands (Ua GLrange style, isotropic)
    size_gl = np.full(len(gl_d), float(COARSE))
    for dist_m, elem_m in GL_BANDS:
        size_gl = np.where(gl_d < dist_m, np.minimum(size_gl, elem_m), size_gl)
    last_dist, last_size = GL_BANDS[-1]
    gl_frac = np.clip((gl_d - last_dist) / last_dist, 0.0, 1.0)
    size_gl = np.where(
        gl_d >= last_dist,
        last_size + (COARSE - last_size) * gl_frac,
        size_gl,
    )

    # 3. Shelf / buffer / ice front
    cf_frac = np.clip(cf_d / 15e3, 0.0, 1.0)
    size_cf_shelf = SHELF_SIZE + (COARSE - SHELF_SIZE) * cf_frac
    cf_frac_buf = np.clip(cf_d / 20e3, 0.0, 1.0)
    size_cf_buf = 8000 + (BUFFER_SIZE - 8000) * cf_frac_buf
    size_cf = np.where(on_shelf, size_cf_shelf,
                       np.where(on_ice, size_cf_shelf, size_cf_buf))

    target_sizes = np.minimum(np.minimum(size_sr, size_gl), size_cf)
    target_sizes = np.clip(target_sizes, GL_BANDS[0][1], COARSE)

    n_1k = int((target_sizes <= 1000).sum())
    n_2k = int((target_sizes <= 2000).sum())
    print(f"  GL bands: <=1km: {n_1k}, <=2km: {n_2k}")
    print(f"  Size range: [{target_sizes.min():.0f}, {target_sizes.max():.0f}] m")

    sf_view = gmsh.view.add("size field")
    gmsh.view.addModelData(sf_view, 0, fn_base + "_raw", "ElementData",
                           tri_tags, target_sizes[:, None])
    gmsh.view.write(sf_view, fn_base + "_sf.pos")

    # Pass 2: remesh with background field + Netgen
    print("\nPass 2: remesh...")
    gmsh.model.add(fn_base)
    build_gmsh_geometry(boundaries, names, GL_BANDS[0][1], COARSE)

    bg = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(bg, "ViewTag", sf_view)
    gmsh.model.mesh.field.setAsBackgroundMesh(bg)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)
    print("  Netgen optimization...")
    gmsh.model.mesh.optimize("Netgen")

    final_tri, _ = gmsh.model.mesh.getElementsByType(2)
    final_verts, _, _ = gmsh.model.mesh.getNodes()
    print(f"  Final: {len(final_verts)} vertices, {len(final_tri)} cells")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(fn_base + ".msh")
    gmsh.finalize()

    for ext in ["_raw.msh", "_sf.pos"]:
        try:
            os.remove(fn_base + ext)
        except FileNotFoundError:
            pass

    size_mb = os.path.getsize(fn_base + ".msh") / 1e6
    print(f"\nSaved: {fn_base}.msh ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
