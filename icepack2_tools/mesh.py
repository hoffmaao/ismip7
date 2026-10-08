#!/usr/bin/env python3
"""
Generate adaptive meshes for Antarctica ice sheet/shelf modeling.

Prerequisites:
    - Run download_data.py first to get BedMachine and velocity data
    - Requires: gmsh, numpy, xarray, rasterio, shapely, geopandas,
                scipy, matplotlib

Workflow:
    1. Extract ice outline from BedMachine mask
    2. Classify boundary segments (calving front vs. other)
    3. Compute velocity-based size field for adaptive refinement
    4. Generate 2D meshes at multiple resolution levels using gmsh

Note: Meshes are 2D (no vertical extrusion) for use with icepack2's
SSA formulation, which handles both grounded and floating ice.

Usage:
    python mesh_antarctica.py
"""

import os
import sys
import glob
import numpy as np
import xarray as xr
import geopandas as gpd
from shapely.geometry import shape, Polygon, MultiPolygon, LineString
from shapely.ops import unary_union
from rasterio.features import shapes as rio_shapes
from rasterio.transform import Affine
from scipy.ndimage import binary_erosion, binary_dilation
import gmsh
import matplotlib.pyplot as plt

# ── Configuration ──────────────────────────────────────────────────────

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# <repo>/data is a gitignored symlink to antarctica/data that a fresh clone does
# not have, so name the observational root the rest of the code uses.
from .runconfig import buffer_m as _buffer_m, obs_data_root
DATA_DIR = obs_data_root()
MESH_DIR = os.path.join(_ROOT, "mesh")

# Outline extraction parameters
SUBSAMPLE = 8  # Subsample BedMachine mask: 500m * 8 = 4 km effective
SIMPLIFY_TOL = 2000.0  # Douglas-Peucker simplification tolerance (m)
MIN_HOLE_AREA = 5e9  # Min nunatak area to keep in mesh (m²), ~70x70 km
MIN_POLY_AREA = 1e10  # Min ice body area to keep (m²), ~100x100 km

# Resolution levels to generate: (fine_target, coarse_target) in meters
# Standard uniform-ratio meshes
FINE_TARGETS = [32000, 16000, 8000, 4000]
# Adaptive mesh: 500m at grounding zone, 32km interior, strain-rate transition
ADAPTIVE_MESH = (500, 32000)


# ── Utility ────────────────────────────────────────────────────────────


def find_file(directory, pattern):
    """Find a file matching a glob pattern in a directory."""
    matches = glob.glob(os.path.join(directory, pattern))
    if not matches:
        raise FileNotFoundError(
            f"No file matching '{pattern}' in {directory}. "
            "Have you run download_data.py?"
        )
    return matches[0]


# ── Step 1: Extract Ice Outline ────────────────────────────────────────


def load_bedmachine_mask(data_dir=DATA_DIR):
    """Load BedMachine mask and coordinates, subsampled for efficiency.

    Mask values: 0=ocean, 1=ice-free land, 2=grounded ice,
                 3=floating ice (shelf), 4=Lake Vostok
    """
    fn = find_file(os.path.join(data_dir, "bedmachine"), "*.nc")
    print(f"Loading BedMachine mask from {fn}...")

    ds = xr.open_dataset(fn)
    mask_full = ds["mask"]

    # Subsample for efficiency
    mask = mask_full.values[::SUBSAMPLE, ::SUBSAMPLE]
    x = ds["x"].values[::SUBSAMPLE]
    y = ds["y"].values[::SUBSAMPLE]
    ds.close()

    print(f"  Subsampled mask shape: {mask.shape}")
    print(f"  x range: [{x[0]/1e3:.0f}, {x[-1]/1e3:.0f}] km")
    print(f"  y range: [{y[0]/1e3:.0f}, {y[-1]/1e3:.0f}] km")
    return mask, x, y


def extract_ice_outline(mask, x, y, buffer_m=None):
    """Extract simplified ice outline polygon from BedMachine mask.

    Returns the main ice sheet polygon (largest connected ice body)
    with simplified boundaries and small nunataks removed, pushed outward by
    ``buffer_m`` metres. None takes runconfig.buffer_m(): ISMIP7_BUFFER_M,
    else the production mesh's 20 km, the value the mesh names assume.
    """
    print("Extracting ice outline...")

    # Binary mask: ice = grounded + floating + Lake Vostok
    ice = ((mask >= 2) & (mask <= 4)).astype(np.uint8)

    # Morphological cleaning: remove small isolated features
    struct = np.ones((3, 3), dtype=bool)
    ice = binary_erosion(ice, struct, iterations=1).astype(np.uint8)
    ice = binary_dilation(ice, struct, iterations=1).astype(np.uint8)

    # Build affine transform from coordinates
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])

    # rasterio expects y decreasing (top to bottom in image)
    if dy > 0:
        # y is increasing: flip array for rasterio convention
        ice_for_rio = ice[::-1, :]
        y_top = y[-1]
        transform = Affine(dx, 0, x[0] - dx / 2, 0, -dy, y_top + dy / 2)
    else:
        ice_for_rio = ice
        transform = Affine(dx, 0, x[0] - dx / 2, 0, dy, y[0] - dy / 2)

    # Vectorize the binary mask to polygons
    print("  Vectorizing mask to polygons...")
    all_polys = []
    for geom, value in rio_shapes(ice_for_rio, transform=transform):
        if value == 1:
            poly = shape(geom)
            if poly.is_valid and poly.area > MIN_POLY_AREA:
                all_polys.append(poly)

    print(f"  Found {len(all_polys)} ice polygon(s) above area threshold")
    if not all_polys:
        raise ValueError("No ice polygons found! Check BedMachine data.")

    # Take the largest polygon (main ice sheet + shelves)
    main_poly = max(all_polys, key=lambda p: p.area)
    print(f"  Main polygon area: {main_poly.area / 1e12:.2f} million km²")

    # Simplify the boundary
    simplified = main_poly.simplify(SIMPLIFY_TOL, preserve_topology=True)

    # Remove ALL interior holes (nunataks). icepack2's dual formulation is
    # well-posed at h=0, so we don't need holes — just set thickness to zero
    # at nunataks and the solver handles it correctly.
    n_holes_removed = len(list(simplified.interiors))
    simplified = Polygon(simplified.exterior)

    # Optionally buffer the outline into the ocean so the calving front
    # is interior to the domain (not at the boundary). This allows
    # dropping the calving_terminus BC entirely — icepack2 handles h=0.
    buffer_dist = _buffer_m() if buffer_m is None else float(buffer_m)
    if buffer_dist > 0:
        simplified = simplified.buffer(buffer_dist, resolution=4)
        # Re-simplify after buffering
        simplified = simplified.simplify(SIMPLIFY_TOL, preserve_topology=True)
        # Remove any holes created by buffering
        if simplified.geom_type == "Polygon":
            simplified = Polygon(simplified.exterior)
        else:
            # MultiPolygon — take largest
            simplified = max(simplified.geoms, key=lambda p: p.area)
            simplified = Polygon(simplified.exterior)
        print(f"  Buffered outline by {buffer_dist/1e3:.0f} km into ocean")

    n_pts = len(simplified.exterior.coords)
    print(f"  Simplified: {n_pts} boundary points")
    print(f"  Removed all {n_holes_removed} interior holes (icepack2 handles h=0)")

    return simplified


# ── Step 2: Classify Boundaries ────────────────────────────────────────


def classify_boundaries(outline, mask, x, y):
    """Classify boundary segments as calving front or other.

    Calving front: ice boundary adjacent to ocean (mask=0)
    Other: ice boundary adjacent to land (mask=1) or ice divide

    Returns lists of coordinate arrays and names for each segment.
    """
    print("Classifying boundary segments...")

    # Ensure exterior is counter-clockwise (shapely convention)
    if not outline.exterior.is_ccw:
        coords = np.array(outline.exterior.coords)[::-1]
    else:
        coords = np.array(outline.exterior.coords)

    n_verts = len(coords)
    n_segs = n_verts - 1  # last vertex = first vertex (closed ring)

    # Compute segment midpoints and outward normals
    seg_dx = coords[1:, 0] - coords[:-1, 0]
    seg_dy = coords[1:, 1] - coords[:-1, 1]
    seg_len = np.sqrt(seg_dx**2 + seg_dy**2)
    seg_len = np.maximum(seg_len, 1.0)  # avoid division by zero

    # Outward normal for CCW polygon: rotate segment vector 90° clockwise
    norm_x = seg_dy / seg_len
    norm_y = -seg_dx / seg_len

    # Sample mask just outside each boundary segment
    mid_x = (coords[:-1, 0] + coords[1:, 0]) / 2.0
    mid_y = (coords[:-1, 1] + coords[1:, 1]) / 2.0
    step = float(abs(x[1] - x[0])) * 2  # 2 pixels outward
    outside_x = mid_x + norm_x * step
    outside_y = mid_y + norm_y * step

    # Create xarray DataArray for interpolation
    mask_da = xr.DataArray(
        mask,
        dims=["y", "x"],
        coords={"x": x.astype(float), "y": y.astype(float)},
    )

    outside_mask = mask_da.interp(
        x=xr.DataArray(outside_x, dims="pts"),
        y=xr.DataArray(outside_y, dims="pts"),
        method="nearest",
    ).values

    # Classify: marine boundary where outside is ocean
    is_marine = outside_mask < 0.5  # ocean = 0

    # Group consecutive segments with same classification
    boundaries = []
    names = []
    start_idx = 0
    current_marine = bool(is_marine[0])

    for i in range(1, n_segs):
        if bool(is_marine[i]) != current_marine:
            # Segment changed type: save the current group
            segment_coords = coords[start_idx : i + 1]
            if current_marine:
                name = f"Calving_{len(boundaries)}"
            else:
                name = f"Other_{len(boundaries)}"
            boundaries.append(segment_coords)
            names.append(name)
            start_idx = i
            current_marine = bool(is_marine[i])

    # Final segment (exclude closing vertex — it duplicates coords[0])
    segment_coords = coords[start_idx:-1]
    if current_marine:
        name = f"Calving_{len(boundaries)}"
    else:
        name = f"Other_{len(boundaries)}"
    boundaries.append(segment_coords)
    names.append(name)

    n_calving = sum(1 for n in names if n.startswith("Calving"))
    n_other = sum(1 for n in names if n.startswith("Other"))
    print(
        f"  {len(boundaries)} boundary segments: "
        f"{n_calving} calving, {n_other} other"
    )

    return boundaries, names


# ── Step 3: Velocity-Based Size Field ──────────────────────────────────


def extract_grounding_line(simplify_tol=5000.0):
    """Extract the grounding line contour from BedMachine as a list of LineStrings.

    The grounding line is where height above flotation = 0.
    Returns simplified contour segments suitable for embedding in gmsh.
    """
    print("Extracting grounding line contour...")
    from shapely.geometry import LineString, MultiLineString
    from shapely.ops import linemerge
    import matplotlib.pyplot as plt

    fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    ds = xr.open_dataset(fn)
    mask = ds["mask"].values[::SUBSAMPLE, ::SUBSAMPLE]
    thick = ds["thickness"].values[::SUBSAMPLE, ::SUBSAMPLE]
    bed = ds["bed"].values[::SUBSAMPLE, ::SUBSAMPLE]
    x = ds["x"].values[::SUBSAMPLE]
    y = ds["y"].values[::SUBSAMPLE]
    ds.close()

    # Height above flotation
    rho_I, rho_W = 917.0, 1024.0
    s_float = bed + (rho_W / rho_I) * np.maximum(-bed, 0)
    s = np.maximum(bed + thick, (1 - rho_I / rho_W) * thick)
    haf = s - s_float

    # Only consider where there's actually ice
    ice = (mask >= 2) & (mask <= 4)
    haf_masked = np.where(ice, haf, np.nan)

    # Extract zero contour
    fig_tmp, ax_tmp = plt.subplots()
    cs = ax_tmp.contour(x, y, haf_masked, levels=[0])
    plt.close(fig_tmp)

    # Convert to shapely LineStrings and simplify
    lines = []
    for path in cs.get_paths():
        verts = path.vertices
        if len(verts) > 2:
            line = LineString(verts).simplify(simplify_tol, preserve_topology=True)
            if line.length > simplify_tol * 5:  # skip tiny fragments
                lines.append(line)

    total_pts = sum(len(l.coords) for l in lines)
    print(
        f"  {len(lines)} grounding line segments, {total_pts} points "
        f"(simplified at {simplify_tol/1e3:.0f}km)"
    )

    return lines


def load_calving_front_field(outline, boundaries, names):
    """Compute a calving front proximity field for mesh refinement.

    Returns a refinement indicator that is high within 50km of any calving
    boundary segment, ensuring the calving retreat experiments have
    sufficient resolution.
    """
    print("Computing calving front proximity field...")
    from scipy.ndimage import distance_transform_edt

    fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    ds = xr.open_dataset(fn)
    x = ds["x"].values[::SUBSAMPLE]
    y = ds["y"].values[::SUBSAMPLE]
    ds.close()
    dx = abs(float(x[1] - x[0]))

    # Rasterize calving boundary segments onto the grid
    calving_mask = np.zeros((len(y), len(x)), dtype=bool)
    for seg_coords, name in zip(boundaries, names):
        if not name.startswith("Calving"):
            continue
        for pt in seg_coords:
            ix = np.argmin(np.abs(x - pt[0]))
            iy = np.argmin(np.abs(y - pt[1]))
            if 0 <= ix < len(x) and 0 <= iy < len(y):
                calving_mask[iy, ix] = True
                # Dilate slightly to catch nearby cells
                for di in range(-2, 3):
                    for dj in range(-2, 3):
                        ii, jj = iy + di, ix + dj
                        if 0 <= ii < len(y) and 0 <= jj < len(x):
                            calving_mask[ii, jj] = True

    # Distance from calving front (in meters)
    dist_cells = distance_transform_edt(~calving_mask)
    dist_m = dist_cells * dx

    # Refinement: full resolution within 50km, decay to coarse beyond
    calving_width = 20e3  # refine within 20km of calving front
    cf_indicator = np.exp(-dist_m / calving_width)

    cf_da = xr.DataArray(
        cf_indicator,
        dims=["y", "x"],
        coords={"x": x.astype(float), "y": y.astype(float)},
    )

    print(f"  Calving boundary pixels: {calving_mask.sum()}")
    print(f"  CF indicator range: [{float(cf_da.min()):.4f}, {float(cf_da.max()):.4f}]")
    return cf_da


def load_grounding_zone_field():
    """Compute a grounding zone proximity field from BedMachine mask.

    Returns a refinement indicator that is high near the grounding line
    (where mask transitions between 2=grounded and 3=floating) and
    decays with distance.
    """
    print("Computing grounding zone proximity field...")
    fn = find_file(os.path.join(DATA_DIR, "bedmachine"), "*.nc")
    ds = xr.open_dataset(fn)
    mask = ds["mask"].values[::SUBSAMPLE, ::SUBSAMPLE]
    x = ds["x"].values[::SUBSAMPLE]
    y = ds["y"].values[::SUBSAMPLE]
    ds.close()

    from scipy.ndimage import distance_transform_edt

    # Grounding zone = boundary between grounded (2) and floating (3)
    grounded = (mask == 2) | (mask == 4)  # grounded + Vostok
    floating = mask == 3

    # Dilate both masks and find overlap = grounding zone region
    struct = np.ones((3, 3), dtype=bool)
    gr_dilated = binary_dilation(grounded, struct, iterations=2)
    fl_dilated = binary_dilation(floating, struct, iterations=2)
    gl_region = gr_dilated & fl_dilated

    # Distance from grounding zone (in grid cells)
    dist_cells = distance_transform_edt(~gl_region)
    dx = abs(float(x[1] - x[0]))
    dist_m = dist_cells * dx

    # Refinement indicator: decays exponentially from grounding zone
    # Lengthscale of ~50km for the transition zone
    gl_indicator = np.exp(-dist_m / 20e3)

    gl_da = xr.DataArray(
        gl_indicator,
        dims=["y", "x"],
        coords={"x": x.astype(float), "y": y.astype(float)},
    )

    print(f"  Grounding zone pixels: {gl_region.sum()}")
    print(f"  GL indicator range: [{float(gl_da.min()):.4f}, {float(gl_da.max()):.4f}]")
    return gl_da


def load_velocity_for_sizing(data_dir=DATA_DIR):
    """Load velocity data and compute a strain-rate-based refinement field.

    The refinement indicator combines strain rate and speed to concentrate
    mesh resolution where ice deformation is high (outlet glaciers, shear
    margins, fast-flowing ice streams).
    """
    print("Loading velocity data for mesh sizing...")

    vel_fn = find_file(os.path.join(data_dir, "velocity"), "*.nc")
    print(f"  Reading {vel_fn}...")

    ds = xr.open_dataset(vel_fn)

    # Handle different variable naming conventions
    vx_key = "VX" if "VX" in ds else "vx"
    vy_key = "VY" if "VY" in ds else "vy"

    # Subsample for efficiency (450m -> ~2km with sub=4)
    sub = 4
    vx = ds[vx_key][::sub, ::sub].astype(np.float64)
    vy = ds[vy_key][::sub, ::sub].astype(np.float64)
    ds.close()

    # Replace invalid values with 0
    vx = vx.where(np.isfinite(vx), 0.0)
    vy = vy.where(np.isfinite(vy), 0.0)

    # Compute approximate strain rate magnitude
    # Using finite differences on the regular grid
    dudx = vx.differentiate("x")
    dudy = vx.differentiate("y")
    dvdx = vy.differentiate("x")
    dvdy = vy.differentiate("y")

    speed = np.sqrt(vx**2 + vy**2)
    mag_eps = np.abs(dudx) + np.abs(dvdy) + 0.5 * (np.abs(dudy) + np.abs(dvdx))

    # Combined refinement indicator (following a collaborator's approach)
    refinement = (mag_eps * 20 + speed / 2000) ** 0.75
    refinement = refinement.fillna(1e-8)
    refinement = refinement.where(refinement > 1e-8, 1e-8)

    print(f"  Velocity grid shape: {vx.shape}")
    print(
        f"  Refinement range: [{float(refinement.min()):.4f}, "
        f"{float(refinement.max()):.4f}]"
    )

    return refinement


# ── Step 4: Mesh Generation ────────────────────────────────────────────


def densify_segment(coords, max_spacing):
    """Add intermediate points so no two consecutive points are more than max_spacing apart."""
    dense = [coords[0]]
    for i in range(1, len(coords)):
        p0 = coords[i - 1]
        p1 = coords[i]
        dist = np.sqrt((p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2)
        if dist > max_spacing:
            n_sub = int(np.ceil(dist / max_spacing))
            for j in range(1, n_sub):
                frac = j / n_sub
                dense.append(p0 + frac * (p1 - p0))
        dense.append(p1)
    return np.array(dense)


def build_gmsh_geometry(boundaries, names, fine_targ, rough_targ):
    """Build gmsh geometry from classified boundary segments.

    Calving front segments get fine_targ resolution at boundary points.
    Calving boundaries are densified to ensure elements are at most fine_targ.
    Other segments get rough_targ resolution at boundary points.
    No interior holes — icepack2 handles h=0 regions natively.
    """
    pt_num = 1
    line_num = 1
    outline_line_groups = []  # list of lists of line tags per segment

    # Add exterior boundary segments
    for seg_idx, seg_coords in enumerate(boundaries):
        is_calving = names[seg_idx].startswith("Calving")
        lc = fine_targ if is_calving else rough_targ

        # Densify calving segments to ensure fine resolution at the boundary
        if is_calving and len(seg_coords) > 1:
            seg_coords = densify_segment(seg_coords, fine_targ)

        seg_lines = []
        first_pt_this_seg = pt_num

        for pt_idx, pt in enumerate(seg_coords):
            # Skip first point of non-first segments (shared with previous)
            if seg_idx > 0 and pt_idx == 0:
                continue
            gmsh.model.geo.addPoint(pt[0], pt[1], 0, lc, pt_num)
            pt_num += 1

        # Build line connectivity
        if seg_idx == 0:
            pts = list(range(first_pt_this_seg, pt_num))
        else:
            # Connect from previous segment's last point
            prev_last_pt = first_pt_this_seg - 1
            pts = [prev_last_pt] + list(range(first_pt_this_seg, pt_num))

        # Close the loop on the last segment
        if seg_idx == len(boundaries) - 1:
            pts.append(1)  # back to first point

        for i in range(len(pts) - 1):
            gmsh.model.geo.addLine(pts[i], pts[i + 1], line_num)
            seg_lines.append(line_num)
            line_num += 1

        outline_line_groups.append(seg_lines)

    # Collect all outline lines for the curve loop
    all_outline_lines = []
    for group in outline_line_groups:
        all_outline_lines.extend(group)

    # Create exterior curve loop
    outline_loop_tag = line_num
    gmsh.model.geo.addCurveLoop(all_outline_lines, outline_loop_tag)
    line_num += 1

    # Create plane surface (no holes — icepack2 handles h=0 natively)
    surface_tag = line_num
    gmsh.model.geo.addPlaneSurface([outline_loop_tag], surface_tag)

    # Physical groups for boundary conditions
    for seg_idx, seg_lines in enumerate(outline_line_groups):
        gmsh.model.geo.addPhysicalGroup(1, seg_lines, name=names[seg_idx])

    gmsh.model.geo.addPhysicalGroup(2, [surface_tag], name="Ice")
    gmsh.model.geo.synchronize()

    return surface_tag, pt_num, line_num


def add_grounding_line_to_gmsh(gl_lines, surface_tag, pt_num, line_num, fine_targ):
    """Embed grounding line contours as interior curves in the gmsh mesh.

    This allows computing flux across the grounding line via dS(gl_id)
    in Firedrake. The GL is densified to fine_targ spacing.

    Returns the physical group tag for the grounding line.
    """
    if not gl_lines:
        return None

    gl_line_tags = []

    for seg_idx, line in enumerate(gl_lines):
        coords = np.array(line.coords)
        # Densify to fine_targ spacing
        coords = densify_segment(coords, fine_targ)
        if len(coords) < 2:
            continue

        seg_pts = []
        for pt in coords:
            gmsh.model.geo.addPoint(pt[0], pt[1], 0, fine_targ, pt_num)
            seg_pts.append(pt_num)
            pt_num += 1

        for i in range(len(seg_pts) - 1):
            gmsh.model.geo.addLine(seg_pts[i], seg_pts[i + 1], line_num)
            gl_line_tags.append(line_num)
            line_num += 1

    if not gl_line_tags:
        return None

    gmsh.model.geo.synchronize()

    # Embed the GL curves in the surface (interior constraint)
    gmsh.model.mesh.embed(1, gl_line_tags, 2, surface_tag)

    # Physical group for the grounding line
    gl_group_tag = gmsh.model.geo.addPhysicalGroup(
        1, gl_line_tags, name="GroundingLine"
    )
    gmsh.model.geo.synchronize()

    print(
        f"    Added grounding line: {len(gl_line_tags)} line segments, "
        f"physical group {gl_group_tag}"
    )

    return gl_group_tag


# ── Marine front of an ice mask (issue #167) ───────────────────────────
#
# A buffered mesh puts the ice front inside the domain. Unless mesh edges
# follow it, the front crosses cells, and vertex sampling gives each of them a
# fraction of the front's thickness: 40 m against BedMachine's 163 m on the
# 2 km mesh (issue #167). The functions below take the marine front of an ice
# mask (obs_icemask: Greene et al. 2022, one year), as curves the mesh
# embeds, so its nodes and edges lie on the front.


def _disk(radius):
    r = int(radius)
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    return (xx * xx + yy * yy) <= r * r


def clean_ice(ice, px, spacing, min_area):
    """The ice a mesh of edge ``spacing`` can follow: channels and spits
    narrower than ``spacing`` closed and opened away, ice islands smaller than
    ``min_area`` [m^2] dropped and holes smaller than it filled."""
    from scipy import ndimage

    disk = _disk(max(1, int(round(spacing / (2.0 * px)))))
    ice = ndimage.binary_closing(np.asarray(ice, bool), structure=disk)
    ice = ndimage.binary_opening(ice, structure=disk)
    npx = float(min_area) / px ** 2
    lab, _ = ndimage.label(ice)
    keep = np.bincount(lab.ravel()) >= npx
    keep[0] = False
    ice = keep[lab]
    lab, _ = ndimage.label(~ice)
    fill = np.bincount(lab.ravel()) < npx
    fill[0] = False
    fill[np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))] = False
    return ice | fill[lab]


def _classes_at(classes, inverse, xy):
    """Nearest-pixel classes at points ``xy``; outside the raster is land."""
    from .obs_icemask import LAND
    col, row = inverse * (xy[:, 0], xy[:, 1])
    col = np.floor(np.asarray(col)).astype(int)
    row = np.floor(np.asarray(row)).astype(int)
    inside = ((row >= 0) & (row < classes.shape[0])
              & (col >= 0) & (col < classes.shape[1]))
    out = np.full(len(xy), LAND, dtype=classes.dtype)
    out[inside] = classes[row[inside], col[inside]]
    return out


def _runs(flags):
    """``(value, start, stop)`` of the runs of a cyclic sequence whose first
    element starts a run."""
    out, start = [], 0
    for i in range(1, len(flags) + 1):
        if i == len(flags) or flags[i] != flags[start]:
            out.append((bool(flags[start]), start, i))
            start = i
    return out


def smooth_marine(flags, lengths, gap, min_run):
    """Edge flags of a closed ring with land runs shorter than ``gap`` turned
    marine and then marine runs shorter than ``min_run`` turned land, both by
    length along the ring."""
    flags = np.asarray(flags, bool).copy()
    lengths = np.asarray(lengths, float)
    for value, short in ((False, gap), (True, min_run)):
        if flags.all() or not flags.any():
            break
        # start at a run boundary, so no run wraps around the end
        shift = int(np.flatnonzero(flags != np.roll(flags, 1))[0])
        f, ln = np.roll(flags, -shift), np.roll(lengths, -shift)
        for v, a, b in _runs(f):
            if v == value and ln[a:b].sum() < short:
                f[a:b] = not value
        flags = np.roll(f, shift)
    return flags


def _resample(line, spacing, closed):
    import shapely
    n = max(3, int(np.ceil(line.length / spacing)))
    t = (np.linspace(0.0, line.length, n, endpoint=False) if closed
         else np.linspace(0.0, line.length, n + 1))
    return shapely.get_coordinates(shapely.line_interpolate_point(line, t))


def _simple(xy, closed):
    from shapely.geometry import LinearRing
    if closed:
        return len(xy) >= 3 and LinearRing(xy).is_simple
    return LineString(xy).is_simple


def extract_marine_front(classes, transform, outline, spacing, min_area=None):
    """The marine front of a classified ice mask as curves for a mesh with
    edge ``spacing`` near the front to embed.

    ``classes`` is :func:`obs_icemask.classify` on a grid with ``transform``.
    The ice is cleaned at the mesh scale (:func:`clean_ice`, ``min_area``
    16 spacing^2 by default) and vectorized; each ring is split, pixel edge by
    pixel edge, into the runs whose ice-free side is marine (land gaps
    shorter than 2 spacing bridged, runs shorter than 3 spacing dropped);
    each run is simplified at one pixel, clipped to 2 spacing inside the mesh
    ``outline``, and resampled at ``spacing``. A curve that is not simple even at half the
    spacing is dropped, and of two curves closer than spacing/2 the shorter
    is.

    Returns ``(curves, stats)``: ``curves`` a list of ``{"xy": (n, 2) array,
    "closed": bool}`` and ``stats`` the counts and lengths.
    """
    import shapely
    from shapely.geometry import LinearRing
    from shapely.geometry.polygon import orient
    from .obs_icemask import ICE, MARINE

    px = abs(float(transform.a))
    min_area = 16.0 * spacing ** 2 if min_area is None else float(min_area)
    ice = clean_ice(classes == ICE, px, spacing, min_area)
    inverse = ~transform
    inner = outline.buffer(-2.0 * spacing)
    stats = dict(spacing=float(spacing), pixel=px, min_area=min_area,
                 polygons=0, rings=0, runs=0, dropped_not_simple=0,
                 dropped_too_close=0, length_km=0.0)

    raw = []
    for geom, value in rio_shapes(ice.astype(np.uint8), mask=ice,
                                  transform=transform, connectivity=4):
        poly = shape(geom)
        if not poly.intersects(inner):
            continue
        # classified on the pixel edges themselves, simplified per run below:
        # a chord of the simplified ring can span several kilometres
        poly = orient(poly, 1.0)
        stats["polygons"] += 1
        for ring in [poly.exterior, *poly.interiors]:
            # one edge per pixel: shapes() merges collinear pixel edges
            xy = densify_segment(np.asarray(ring.coords), px)[:-1]
            if len(xy) < 3:
                continue
            stats["rings"] += 1
            d = np.roll(xy, -1, axis=0) - xy
            length = np.hypot(d[:, 0], d[:, 1])
            normal = np.zeros_like(d)
            ok = length > 0
            # ice on the left of every ring after orient(), so the right
            # normal points at the ice-free side
            normal[ok] = np.column_stack([d[ok, 1], -d[ok, 0]]) / length[ok, None]
            mid = xy + 0.5 * d
            marine = np.zeros(len(xy), bool)
            for k in (1.0, 2.0):
                marine |= _classes_at(classes, inverse, mid + k * px * normal) == MARINE
            marine = smooth_marine(marine, length, 2.0 * spacing, 3.0 * spacing)
            if marine.all():
                raw.append((LinearRing(xy).simplify(px, preserve_topology=True), True))
                continue
            shift = int(np.flatnonzero(marine != np.roll(marine, 1))[0])
            xy_r, m_r = np.roll(xy, -shift, axis=0), np.roll(marine, -shift)
            for v, a, b in _runs(m_r):
                if v:
                    pts = np.vstack([xy_r[a:b], xy_r[b % len(xy_r)]])
                    raw.append((LineString(pts).simplify(px, preserve_topology=True), False))

    pieces = []
    for line, closed in raw:
        if closed and line.within(inner):
            pieces.append((line, True))
            continue
        clipped = LineString(line.coords).intersection(inner)
        for part in getattr(clipped, "geoms", [clipped]):
            if part.geom_type == "LineString" and part.length >= 3.0 * spacing:
                pieces.append((part, False))

    resampled = []
    for line, closed in pieces:
        for h in (spacing, 0.5 * spacing):
            xy = _resample(line, h, closed)
            if _simple(xy, closed):
                resampled.append((xy, closed))
                break
        else:
            stats["dropped_not_simple"] += 1

    geoms = [LinearRing(xy) if closed else LineString(xy) for xy, closed in resampled]
    order = np.argsort([-g.length for g in geoms])
    kept = []
    if geoms:
        tree = shapely.STRtree(geoms)
        near = tree.query(geoms, predicate="dwithin", distance=0.5 * spacing)
        partners = {}
        for i, j in zip(*near):
            if i != j:
                partners.setdefault(int(i), set()).add(int(j))
        chosen = set()
        for i in order:
            if partners.get(int(i), set()) & chosen:
                stats["dropped_too_close"] += 1
                continue
            chosen.add(int(i))
        kept = [resampled[i] for i in sorted(chosen)]
    curves = [{"xy": xy, "closed": closed} for xy, closed in kept]
    stats["runs"] = len(curves)
    stats["closed"] = sum(c["closed"] for c in curves)
    stats["points"] = int(sum(len(c["xy"]) for c in curves))
    stats["length_km"] = float(sum(
        (LinearRing(c["xy"]) if c["closed"] else LineString(c["xy"])).length
        for c in curves)) / 1e3
    return curves, stats


def densify_curves(curves, spacing):
    """Every curve's points with points added so none is more than
    ``spacing`` from the next, all stacked: the points a distance to the
    front is measured from."""
    out = []
    for c in curves:
        xy = c["xy"]
        if c["closed"]:
            xy = np.vstack([xy, xy[:1]])
        out.append(densify_segment(xy, spacing))
    return np.vstack(out) if out else np.zeros((0, 2))


def embed_front_in_gmsh(curves, surface_tag, pt_num, line_num, size):
    """Embed ``curves`` (from :func:`extract_marine_front`) in the surface as
    interior constraints, with no physical group: the mesh's nodes and edges
    then follow the front, and the boundary-id sidecar is unchanged (an
    interior physical line would become a facet marker that
    ``boundary.load_boundary_ids`` refuses). Every curve point becomes a node.

    Returns ``(pt_num, line_num, points)``, the next free tags and the number
    of points added.
    """
    tags, points = [], 0
    line_num += 1                 # build_gmsh_geometry's last tag is in use
    for c in curves:
        first = pt_num
        for x, y in c["xy"]:
            gmsh.model.geo.addPoint(float(x), float(y), 0.0, size, pt_num)
            pt_num += 1
        ids = list(range(first, pt_num))
        if c["closed"]:
            ids.append(first)
        for a, b in zip(ids[:-1], ids[1:]):
            gmsh.model.geo.addLine(a, b, line_num)
            tags.append(line_num)
            line_num += 1
        points += len(c["xy"])
    gmsh.model.geo.synchronize()
    if tags:
        gmsh.model.mesh.embed(1, tags, 2, surface_tag)
    return pt_num, line_num, points


def front_points_off_mesh(curves, nodes_xy, tol=1e-3):
    """How many curve points are not mesh nodes (within ``tol`` metres)."""
    from scipy.spatial import cKDTree
    if not curves:
        return 0
    pts = np.vstack([c["xy"] for c in curves])
    d, _ = cKDTree(np.asarray(nodes_xy)[:, :2]).query(pts)
    return int((d > tol).sum())


def save_front(curves, path, stats=None):
    """Write curves to an .npz (points stacked, with each curve's start and
    closed flag), the record of what a mesh was built to follow."""
    xy = np.vstack([c["xy"] for c in curves]) if curves else np.zeros((0, 2))
    starts = np.cumsum([0] + [len(c["xy"]) for c in curves])
    np.savez(path, xy=xy, starts=starts,
             closed=np.array([c["closed"] for c in curves], bool),
             stats=np.array(repr(stats or {})))


def load_front(path):
    """Curves written by :func:`save_front`."""
    with np.load(path) as f:
        xy, starts, closed = f["xy"], f["starts"], f["closed"]
    return [{"xy": xy[a:b], "closed": bool(c)}
            for a, b, c in zip(starts[:-1], starts[1:], closed)]


def generate_mesh(fine_targ, boundaries, names, refinement):
    """Generate a mesh at the given resolution using gmsh's two-pass approach.

    Pass 1: Generate raw mesh with boundary-prescribed sizes.
    Pass 2: Re-mesh with a background size field computed from velocity
             strain rates on the raw mesh.
    """
    rough_targ = fine_targ * 10
    fn_base = os.path.join(MESH_DIR, f"antarctica_{rough_targ}_{fine_targ}")

    print(
        f"\nGenerating mesh: fine={fine_targ/1e3:.0f} km, "
        f"coarse={rough_targ/1e3:.0f} km"
    )

    # ── Pass 1: Raw mesh ──
    gmsh.initialize(sys.argv)
    gmsh.option.setNumber("General.Verbosity", 2)  # reduce output noise
    gmsh.model.add(fn_base + "_raw")

    build_gmsh_geometry(boundaries, names, fine_targ, rough_targ)
    gmsh.model.mesh.generate(2)
    gmsh.write(fn_base + "_raw.msh")

    # Extract mesh data for size field computation
    vtags, vxyz, _ = gmsh.model.mesh.getNodes()
    vxyz = vxyz.reshape((-1, 3))
    vmap = {int(j): i for i, j in enumerate(vtags)}

    tri_tags, tri_vtags = gmsh.model.mesh.getElementsByType(2)
    tri_vids = np.array([vmap[int(j)] for j in tri_vtags])
    triangles = tri_vids.reshape((tri_tags.shape[-1], -1))

    # Compute element-wise target sizes from refinement field
    tri_centers = vxyz[triangles].mean(axis=1)

    ref_x = refinement.coords[refinement.dims[-1]]  # x coordinate
    ref_y = refinement.coords[refinement.dims[-2]]  # y coordinate

    cx = xr.DataArray(tri_centers[:, 0], dims="tri")
    cy = xr.DataArray(tri_centers[:, 1], dims="tri")

    ref_vals = refinement.interp(
        {ref_x.name: cx, ref_y.name: cy}, method="nearest"
    ).values.flatten()
    ref_vals = np.nan_to_num(ref_vals, nan=1e-8)
    ref_vals = np.maximum(ref_vals, 1e-8)

    target_sizes = fine_targ / ref_vals
    target_sizes = np.clip(target_sizes, fine_targ, rough_targ)

    # Create size field view
    sf_view = gmsh.view.add("size field")
    gmsh.view.addModelData(
        sf_view, 0, fn_base + "_raw", "ElementData", tri_tags, target_sizes[:, None]
    )
    gmsh.view.write(sf_view, fn_base + "_sf.pos")

    print(f"  Pass 1 done: {len(tri_tags)} elements")

    # ── Pass 2: Re-mesh with background size field ──
    gmsh.model.add(fn_base)

    build_gmsh_geometry(boundaries, names, fine_targ, rough_targ)

    # Apply background mesh size field
    bg_field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)

    # Get final element count
    final_tri_tags, _ = gmsh.model.mesh.getElementsByType(2)
    print(f"  Pass 2 done: {len(final_tri_tags)} elements")

    # Save as gmsh format 2.2 (compatible with Firedrake)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(fn_base + ".msh")

    gmsh.finalize()

    # Clean up intermediate files
    for ext in ["_raw.msh", "_sf.pos"]:
        try:
            os.remove(fn_base + ext)
        except FileNotFoundError:
            pass

    print(f"  Saved: {fn_base}.msh")


def generate_adaptive_mesh(
    fine_targ,
    coarse_targ,
    boundaries,
    names,
    refinement,
    gl_field,
    cf_field=None,
    gl_lines=None,
):
    """Generate adaptive mesh with grounding-zone + strain-rate refinement.

    Combines grounding zone proximity (500m resolution) with strain-rate
    based refinement, transitioning smoothly to coarse_targ in the interior.
    """
    fn_base = os.path.join(MESH_DIR, f"antarctica_{coarse_targ}_{fine_targ}")

    print(
        f"\nGenerating adaptive mesh: fine={fine_targ}m, "
        f"coarse={coarse_targ/1e3:.0f}km"
    )

    # ── Pass 1: Raw mesh at intermediate resolution ──
    intermediate = min(coarse_targ, 8000)
    gmsh.initialize(sys.argv)
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add(fn_base + "_raw")

    build_gmsh_geometry(boundaries, names, intermediate, coarse_targ)
    gmsh.model.mesh.generate(2)
    gmsh.write(fn_base + "_raw.msh")

    vtags, vxyz, _ = gmsh.model.mesh.getNodes()
    vxyz = vxyz.reshape((-1, 3))
    vmap = {int(j): i for i, j in enumerate(vtags)}

    tri_tags, tri_vtags = gmsh.model.mesh.getElementsByType(2)
    tri_vids = np.array([vmap[int(j)] for j in tri_vtags])
    triangles = tri_vids.reshape((tri_tags.shape[-1], -1))

    tri_centers = vxyz[triangles].mean(axis=1)

    # Interpolate strain-rate refinement
    ref_x = refinement.coords[refinement.dims[-1]]
    ref_y = refinement.coords[refinement.dims[-2]]
    cx = xr.DataArray(tri_centers[:, 0], dims="tri")
    cy = xr.DataArray(tri_centers[:, 1], dims="tri")

    ref_vals = refinement.interp(
        {ref_x.name: cx, ref_y.name: cy}, method="nearest"
    ).values.flatten()
    ref_vals = np.nan_to_num(ref_vals, nan=1e-8)
    ref_vals = np.maximum(ref_vals, 1e-8)

    # Interpolate grounding zone proximity
    gl_x = gl_field.coords[gl_field.dims[-1]]
    gl_y = gl_field.coords[gl_field.dims[-2]]
    gl_vals = gl_field.interp(
        {gl_x.name: cx, gl_y.name: cy}, method="nearest"
    ).values.flatten()
    gl_vals = np.nan_to_num(gl_vals, nan=0.0)

    # Combined size field:
    # - Strain-rate: refines to max(fine_targ*4, 8km) in high-deformation areas
    #   (NOT all the way to fine_targ — that's reserved for GL/calving proximity)
    # - GL proximity: refines to fine_targ near grounding zone
    # - CF proximity: refines to fine_targ near calving front
    strainrate_floor = max(fine_targ * 4, 8000)  # strain rate doesn't go below 8km
    size_strainrate = np.clip(
        strainrate_floor / ref_vals, strainrate_floor, coarse_targ
    )
    size_gl = fine_targ + (coarse_targ - fine_targ) * (1.0 - gl_vals)
    target_sizes = np.minimum(size_strainrate, size_gl)

    if cf_field is not None:
        cf_x = cf_field.coords[cf_field.dims[-1]]
        cf_y = cf_field.coords[cf_field.dims[-2]]
        cf_vals = cf_field.interp(
            {cf_x.name: cx, cf_y.name: cy}, method="nearest"
        ).values.flatten()
        cf_vals = np.nan_to_num(cf_vals, nan=0.0)
        size_cf = fine_targ + (coarse_targ - fine_targ) * (1.0 - cf_vals)
        target_sizes = np.minimum(target_sizes, size_cf)

    target_sizes = np.clip(target_sizes, fine_targ, coarse_targ)

    # Create size field view
    sf_view = gmsh.view.add("adaptive size field")
    gmsh.view.addModelData(
        sf_view, 0, fn_base + "_raw", "ElementData", tri_tags, target_sizes[:, None]
    )
    gmsh.view.write(sf_view, fn_base + "_sf.pos")

    print(f"  Pass 1 done: {len(tri_tags)} elements")
    print(
        f"  Target size range: [{target_sizes.min():.0f}, {target_sizes.max():.0f}] m"
    )

    # ── Pass 2: Re-mesh with background size field + grounding line ──
    gmsh.model.add(fn_base)
    surface_tag, pt_num, line_num = build_gmsh_geometry(
        boundaries, names, fine_targ, coarse_targ
    )

    # Embed grounding line as interior curve (for flux computation)
    gl_group = None
    if gl_lines is not None:
        gl_group = add_grounding_line_to_gmsh(
            gl_lines, surface_tag, pt_num, line_num, fine_targ
        )

    bg_field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)

    final_tri_tags, _ = gmsh.model.mesh.getElementsByType(2)
    print(f"  Pass 2 done: {len(final_tri_tags)} elements")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(fn_base + ".msh")

    gmsh.finalize()

    for ext in ["_raw.msh", "_sf.pos"]:
        try:
            os.remove(fn_base + ext)
        except FileNotFoundError:
            pass

    print(f"  Saved: {fn_base}.msh")
    if gl_group is not None:
        print(f"  Grounding line physical group: {gl_group}")


# ── Visualization ──────────────────────────────────────────────────────


def plot_outline(outline, boundaries, names, save_path=None):
    """Plot the extracted outline with classified boundaries."""
    fig, ax = plt.subplots(figsize=(10, 10))

    for seg_coords, name in zip(boundaries, names):
        color = "red" if name.startswith("Calving") else "black"
        lw = 1.5 if name.startswith("Calving") else 0.5
        ax.plot(
            seg_coords[:, 0] / 1e3,
            seg_coords[:, 1] / 1e3,
            color=color,
            linewidth=lw,
        )

    ax.set_aspect("equal")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_title(
        "Antarctica Ice Outline (EPSG:3031)\n"
        "Red = calving front, Black = other boundary"
    )

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"  Saved outline plot: {save_path}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────────────────


def main():
    os.makedirs(MESH_DIR, exist_ok=True)
    os.makedirs(os.path.join(_ROOT, "figs"), exist_ok=True)

    # Step 1: Extract outline from BedMachine
    mask, x, y = load_bedmachine_mask()
    # The meshes below carry no _buffered tag, the naming of the unbuffered
    # legacy meshes; mesh_antarctica.py builds the named, buffered ones.
    outline = extract_ice_outline(mask, x, y, buffer_m=0.0)

    # Save outline for inspection/reuse
    gdf = gpd.GeoDataFrame(geometry=[outline], crs="EPSG:3031")
    outline_fn = os.path.join(MESH_DIR, "antarctica_ice_outline.gpkg")
    gdf.to_file(outline_fn, driver="GPKG")
    print(f"  Saved outline: {outline_fn}")

    # Step 2: Classify boundary segments
    boundaries, names = classify_boundaries(outline, mask, x, y)

    # Plot outline
    plot_outline(outline, boundaries, names, save_path="figs/antarctica_outline.png")

    # Step 3: Load velocity + grounding zone + calving front for mesh sizing
    refinement = load_velocity_for_sizing()
    gl_field = load_grounding_zone_field()
    cf_field = load_calving_front_field(outline, boundaries, names)

    # Step 4: Generate meshes
    # Standard uniform-ratio meshes
    for fine_targ in FINE_TARGETS:
        generate_mesh(fine_targ, boundaries, names, refinement)

    # Extract grounding line contour for interior mesh constraint
    gl_lines = extract_grounding_line(simplify_tol=2000.0)

    # Adaptive mesh: 500m grounding zone + calving front, 32km interior
    # with grounding line as interior boundary for flux computation
    fine, coarse = ADAPTIVE_MESH
    generate_adaptive_mesh(
        fine, coarse, boundaries, names, refinement, gl_field, cf_field, gl_lines
    )

    print("\n" + "=" * 60)
    print("All meshes generated!")
    for fn in sorted(glob.glob(os.path.join(MESH_DIR, "*.msh"))):
        size_mb = os.path.getsize(fn) / 1e6
        print(f"  {fn}  ({size_mb:.1f} MB)")
    print("=" * 60)


if __name__ == "__main__":
    main()
