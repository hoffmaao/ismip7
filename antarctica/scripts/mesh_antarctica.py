#!/usr/bin/env python3
"""
Generate anisotropic Antarctic mesh with grounding zone refinement.

Uses gmsh BAMG with a metric tensor field for anisotropic elements
at the grounding line: fine perpendicular to the GL (resolves the
grounded/floating transition), coarse along it.

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

# Anisotropic targets
GL_PERP = 500  # 500m perpendicular to GL
GL_PAR = 15000  # 15km along GL
GL_RANGE = 3e3  # anisotropy within 3km of GL
COARSE = 64000  # 64km interior


def find_file(d, p):
    m = glob.glob(os.path.join(d, p))
    if not m:
        raise FileNotFoundError(f"No {p} in {d}")
    return m[0]


def compute_gl_distance_and_normal():
    """Compute GL distance and normal direction from BedMachine."""
    from scipy.ndimage import distance_transform_edt, binary_dilation

    print("Computing GL distance + normal...")
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
    gr_dil = binary_dilation(grounded, struct, iterations=2)
    fl_dil = binary_dilation(floating, struct, iterations=2)
    gl_region = gr_dil & fl_dil & ice

    dist = distance_transform_edt(~gl_region) * dx

    # GL normal: gradient of HAF (points from floating → grounded)
    # Smooth HAF first
    from scipy.ndimage import gaussian_filter

    haf_smooth = gaussian_filter(haf.astype(float), sigma=3)
    grad_y, grad_x = np.gradient(haf_smooth, dx)
    mag = np.sqrt(grad_x**2 + grad_y**2) + 1e-10
    nx = grad_x / mag  # normal x component
    ny = grad_y / mag  # normal y component

    has_ice = (thick > 10) & ice
    is_float = mask_bm == 3

    dist_da = xr.DataArray(
        dist, dims=["y", "x"], coords={"x": x.astype(float), "y": y.astype(float)}
    )
    nx_da = xr.DataArray(
        nx, dims=["y", "x"], coords={"x": x.astype(float), "y": y.astype(float)}
    )
    ny_da = xr.DataArray(
        ny, dims=["y", "x"], coords={"x": x.astype(float), "y": y.astype(float)}
    )
    ice_da = xr.DataArray(
        has_ice.astype(float),
        dims=["y", "x"],
        coords={"x": x.astype(float), "y": y.astype(float)},
    )
    float_da = xr.DataArray(
        is_float.astype(float),
        dims=["y", "x"],
        coords={"x": x.astype(float), "y": y.astype(float)},
    )

    print(f"  GL pixels: {gl_region.sum()}")
    return dist_da, nx_da, ny_da, ice_da, float_da


def build_metric_tensor(h_perp, h_par, nx, ny):
    """Build 2D metric tensor from perpendicular/parallel sizes and normal.

    M = (1/h_perp²) n⊗n + (1/h_par²) t⊗t
    where t = (-ny, nx) is the tangent.

    Returns m11, m12, m22 (symmetric 2D tensor components).
    """
    # n⊗n
    nn11 = nx * nx
    nn12 = nx * ny
    nn22 = ny * ny

    # t⊗t (t = (-ny, nx))
    tt11 = ny * ny
    tt12 = -nx * ny
    tt22 = nx * nx

    m11 = nn11 / h_perp**2 + tt11 / h_par**2
    m12 = nn12 / h_perp**2 + tt12 / h_par**2
    m22 = nn22 / h_perp**2 + tt22 / h_par**2

    return m11, m12, m22


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

    gl_dist, nx_field, ny_field, ice_field, float_field = (
        compute_gl_distance_and_normal()
    )

    fn_base = os.path.join(MESH_DIR, f"antarctica_{COARSE}_{GL_PERP}_aniso")

    # ── Pass 1: raw isotropic mesh ──
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
    n_tri = len(tri_tags)

    print(f"  {len(vtags)} nodes, {n_tri} triangles")

    # ── Interpolate fields at triangle vertices ──
    # (BAMG needs nodal metric, not element-center)
    print("  Interpolating fields at vertices...")
    vx = xr.DataArray(vxyz[:, 0], dims="v")
    vy = xr.DataArray(vxyz[:, 1], dims="v")

    def interp_field(da):
        fx = da.coords[da.dims[-1]]
        fy = da.coords[da.dims[-2]]
        return np.nan_to_num(
            da.interp({fx.name: vx, fy.name: vy}, method="nearest").values.flatten(),
            nan=0.0,
        )

    gl_d = interp_field(gl_dist)
    gl_nx = interp_field(nx_field)
    gl_ny = interp_field(ny_field)
    on_ice = interp_field(ice_field) > 0.5
    on_shelf = interp_field(float_field) > 0.5

    ref_x = refinement.coords[refinement.dims[-1]]
    ref_y = refinement.coords[refinement.dims[-2]]
    ref_vals = np.maximum(
        np.nan_to_num(
            refinement.interp(
                {ref_x.name: vx, ref_y.name: vy}, method="nearest"
            ).values.flatten(),
            nan=1e-8,
        ),
        1e-8,
    )

    # ── Build metric at each vertex ──
    print("  Building metric tensor field...")

    # Isotropic base: strain-rate refinement + moderate shelf/buffer sizing
    sr_floor = 8000
    h_iso = np.clip(sr_floor / ref_vals, sr_floor, COARSE)
    # Shelves get 2km at the front, buffer gets 10km
    h_iso = np.where(on_shelf, np.minimum(h_iso, 5000), h_iso)
    h_iso = np.where(~on_ice, np.minimum(h_iso, 10000), h_iso)

    # Anisotropic GL refinement: ONLY on grounded ice (not buffer/shelf)
    on_grounded = on_ice & ~on_shelf
    gl_frac = np.clip(gl_d / GL_RANGE, 0.0, 1.0)  # 0 at GL, 1 at GL_RANGE

    # Perpendicular: GL_PERP near GL, ramp to h_iso far away
    # Only apply on grounded ice; elsewhere stay isotropic
    h_perp = np.where(on_grounded, GL_PERP + (h_iso - GL_PERP) * gl_frac, h_iso)
    # Parallel: GL_PAR near GL, ramp to h_iso far away
    h_par = np.where(on_grounded, GL_PAR + (h_iso - GL_PAR) * gl_frac, h_iso)

    # Ensure h_perp <= h_par (anisotropy only stretches, doesn't compress)
    h_perp = np.minimum(h_perp, h_par)
    # Ensure minimum sizes
    h_perp = np.maximum(h_perp, GL_PERP)
    h_par = np.maximum(h_par, GL_PERP)

    # Build metric: M = (1/hp²) n⊗n + (1/hq²) t⊗t
    m11, m12, m22 = build_metric_tensor(h_perp, h_par, gl_nx, gl_ny)

    n_aniso = int((h_par / h_perp > 2).sum())
    print(f"  Anisotropic nodes (ratio > 2): {n_aniso}")
    print(f"  h_perp range: [{h_perp.min():.0f}, {h_perp.max():.0f}] m")
    print(f"  h_par range:  [{h_par.min():.0f}, {h_par.max():.0f}] m")

    # ── Write metric as TT .pos file ──
    print("  Writing metric .pos file...")
    pos_fn = fn_base + "_metric.pos"
    with open(pos_fn, "w") as f:
        f.write('View "nodalMetric" {\n')
        for t_idx in range(n_tri):
            v0, v1, v2 = triangles[t_idx]
            # Coordinates
            coords_str = ",".join(
                f"{vxyz[v, 0]},{vxyz[v, 1]},{vxyz[v, 2]}" for v in [v0, v1, v2]
            )
            # Tensor values: m11,m12,0, m12,m22,0, 0,0,1 per node
            tensors = []
            for v in [v0, v1, v2]:
                tensors.extend([m11[v], m12[v], 0, m12[v], m22[v], 0, 0, 0, 1])
            tensor_str = ",".join(f"{t:.8e}" for t in tensors)
            f.write(f"TT({coords_str}){{{tensor_str}}};\n")
        f.write("};\n")

    pos_size = os.path.getsize(pos_fn) / 1e6
    print(f"  Metric file: {pos_fn} ({pos_size:.1f} MB)")

    # ── Pass 2: anisotropic re-mesh with BAMG ──
    print("\nPass 2: BAMG anisotropic mesh...")
    gmsh.model.add(fn_base)
    build_gmsh_geometry(boundaries, names, 8000, COARSE)  # coarse boundary

    # Load metric
    gmsh.merge(pos_fn)
    bg = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(bg, "ViewIndex", 0)
    gmsh.model.mesh.field.setAsBackgroundMesh(bg)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 7)  # BAMG
    gmsh.option.setNumber("Mesh.AnisoMax", 50)  # max anisotropy ratio
    gmsh.option.setNumber("Mesh.SmoothRatio", 3)

    gmsh.model.mesh.generate(2)

    # Netgen optimization
    print("  Netgen optimization...")
    gmsh.model.mesh.optimize("Netgen")

    final_tri, _ = gmsh.model.mesh.getElementsByType(2)
    final_verts, _, _ = gmsh.model.mesh.getNodes()
    print(f"  Final: {len(final_verts)} vertices, {len(final_tri)} cells")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(fn_base + ".msh")
    gmsh.finalize()

    # Clean up
    for ext in ["_raw.msh", "_metric.pos"]:
        try:
            os.remove(fn_base + ext)
        except FileNotFoundError:
            pass

    size_mb = os.path.getsize(fn_base + ".msh") / 1e6
    print(f"\nSaved: {fn_base}.msh ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
