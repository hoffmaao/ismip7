#! /usr/bin/env python3
# vim:fenc=utf-8
#
# Copyright © 2026 David Lilien <dlilien@iu.edu>
#
"""
Build Greenland meshes from a velocity mosaic and an outline.

Velocity mosaics are used for mesh size and outline determines extent.

Usage:
    mesh_greenland.py --data-dir <data_dir> --outline-file <outline_file> --stride <stride> --output-dir <output_dir> --min-lc <min_lc> --num-levels <num_levels> --plot-level <plot_level> --no-plot --boundaries-file <boundaries_file>

    mesh_greenland.py --help
"""
import argparse
from pathlib import Path
from typing import Literal, Sequence

import geopandas as gpd
import gmsh
import numpy as np
from numpy.typing import NDArray
import rioxarray as rxr
import xarray as xr
from shapely import simplify
from shapely.geometry import LineString, MultiPolygon, Polygon


EPSG_GREENLAND = "EPSG:3413"
VX_FILENAME = "greenland_vel_mosaic250_vx_v1.tif"
VY_FILENAME = "greenland_vel_mosaic250_vy_v1.tif"
PROMICE_OUTLINE_FILENAME = "02-PROMICE-2022-IceMask-polygon-v3.gpkg"
SIMPLE_OUTLINE_FILENAME = "simple_polygon_of_greenland.gpkg"
BUFFERED_OUTLINE_FILENAME = "greenland_buffered_5km_simp.gpkg"
DETAILED_OUTLINE_FILENAME = "greenland_outline_mod_from_IMBIE.gpkg"
MULTIYEAR_VELOCITY_DIR = Path("velocity") / "multiyear"
PROMICE_OUTLINE_DIR = "promice"
CUSTOM_BOUNDARIES_DIR = "custom_boundaries"
PROMICE_SIMPLIFY_TOLERANCE = 250.0
DEFAULT_FIG_DIR = Path(__file__).resolve().parent.parent / "figs"

OutlineKind = Literal["promice", "simple", "detailed", "buffered"]
OUTLINE_KINDS: tuple[OutlineKind, ...] = ("promice", "simple", "detailed", "buffered")
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.integer]


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Build Greenland meshes from a velocity mosaic and an outline. "
            "Velocity mosaics are used for mesh size and outline determines extent."
        )
    )
    parser.add_argument(
        "--data-dir",
        "--data_dir",
        dest="data_dir",
        type=Path,
        help=(
            "Directory containing the download_data.py layout. The multiyear "
            "velocity mosaic is read from velocity/multiyear, and the PROMICE "
            "geopackage is read from promice if --outline-file is omitted."
        ),
        default=Path("../data/"),
    )
    parser.add_argument(
        "--outline-kind",
        choices=OUTLINE_KINDS,
        default="promice",
        help=(
            "Built-in outline to mesh. Choices map to the PROMICE geopackage, "
            "custom_boundaries/simple_polygon_of_greenland.gpkg, "
            "custom_boundaries/greenland_buffered_5km_simp.gpkg, and "
            "custom_boundaries/greenland_outline_mod_from_IMBIE.gpkg."
        ),
    )
    parser.add_argument(
        "--outline-file",
        type=Path,
        default=None,
        help="Custom outline geopackage to mesh. Use with --outline-name.",
    )
    parser.add_argument(
        "--outline-name",
        default=None,
        help="Custom outline name used in output filenames. Use with --outline-file.",
    )
    parser.add_argument(
        "--mesh-kind",
        choices=OUTLINE_KINDS,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=None,
        help=(
            "Number of segments to split each outline edge into. If omitted, "
            "the default is inferred from the outline filename: PROMICE=1, "
            "simple=32, otherwise=2."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("../meshes/"),
        help="Directory for generated meshes. Defaults to ../meshes/",
    )
    parser.add_argument(
        "--outline-simplify",
        type=float,
        default=None,
        help=(
            "Simplify the outline by this tolerance before meshing. Defaults "
            "to 250 m for PROMICE outlines and no simplification otherwise."
        ),
    )
    parser.add_argument(
        "--cutoff-velocity",
        type=float,
        default=10.0,
        help="Normal velocity threshold used to tag outlet boundaries.",
    )
    parser.add_argument(
        "--min-lc",
        type=int,
        default=250,
        help="Smallest fine target length scale.",
    )
    parser.add_argument(
        "--num-levels",
        type=int,
        default=8,
        help="Number of mesh resolution levels to build.",
    )
    parser.add_argument(
        "--plot-level",
        nargs=2,
        type=int,
        default=(20000, 2000),
        metavar=("ROUGH", "FINE"),
        help="Mesh level to plot after generation.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip writing the final Firedrake triplot PNG.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=DEFAULT_FIG_DIR,
        help="Directory for mesh figures.",
    )
    parser.add_argument(
        "--boundaries-file",
        type=Path,
        default=Path("boundaries_greenland.gpkg"),
        help="Geopackage path for outlet/other boundary linework.",
    )
    return parser.parse_args()


def candidate_outline_files(data_dir: Path, outline_kind: OutlineKind) -> tuple[Path, ...]:
    """
    Return candidate files for a built-in outline kind.

    Parameters
    ----------
    data_dir : Path
        Path to the data directory.
    outline_kind : {"promice", "simple", "detailed", "buffered"}
        Built-in outline kind.

    Returns
    -------
    tuple[pathlib.Path, ...]
        Candidate outline files in priority order.
    """
    if outline_kind == "promice":
        return (
            data_dir / PROMICE_OUTLINE_DIR / PROMICE_OUTLINE_FILENAME,
            data_dir / PROMICE_OUTLINE_FILENAME,
        )
    filenames = {
        "simple": SIMPLE_OUTLINE_FILENAME,
        "buffered": BUFFERED_OUTLINE_FILENAME,
        "detailed": DETAILED_OUTLINE_FILENAME,
    }
    filename = filenames[outline_kind]
    return (
        data_dir / CUSTOM_BOUNDARIES_DIR / filename,
        data_dir / filename,
    )


def outline_file_for_kind(data_dir: Path, outline_kind: OutlineKind) -> Path:
    """
    Determine the outline file for a built-in outline kind.

    Parameters
    ----------
    data_dir : Path
        Path to the data directory.
    outline_kind : {"promice", "simple", "detailed", "buffered"}
        Built-in outline kind.

    Returns
    -------
    pathlib.Path
        Path to the selected outline file.
    """
    candidates = candidate_outline_files(data_dir, outline_kind)
    for outline_file in candidates:
        if outline_file.exists():
            return outline_file
    candidates_text = " or ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"No {outline_kind} outline file found at {candidates_text}")


def validate_outline_name(outline_name: str) -> str:
    """
    Validate an outline name for use in output filenames.

    Parameters
    ----------
    outline_name : str
        Outline name to validate.

    Returns
    -------
    str
        Validated outline name.
    """
    if not outline_name or "/" in outline_name or "\\" in outline_name:
        raise ValueError("outline name must be nonempty and cannot contain slashes")
    return outline_name


def resolve_outline(args: argparse.Namespace, data_dir: Path) -> tuple[Path, str]:
    """
    Resolve the outline file and output name from command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    data_dir : Path
        Path to the data directory.

    Returns
    -------
    tuple[pathlib.Path, str]
        Outline file and output name.
    """
    if args.outline_file is None and args.outline_name is None:
        outline_kind_name = args.mesh_kind if args.mesh_kind is not None else args.outline_kind
        return outline_file_for_kind(data_dir, outline_kind_name), outline_kind_name

    if args.outline_file is None or (args.outline_name is None and args.mesh_kind is None):
        raise ValueError(
            "Use --outline-file together with --outline-name for custom outlines"
        )

    outline_name = args.outline_name if args.outline_name is not None else args.mesh_kind
    return args.outline_file.expanduser(), validate_outline_name(outline_name)


def outline_kind(outline_file: Path) -> OutlineKind:
    """
    Determine naming convention for the outline file.

    Also used to determine the default stride.

    Parameters
    ----------
    outline_file : Path
        Path to the outline file.

    Returns
    -------
    {"promice", "simple", "detailed", "buffered"}
        Outline type inferred from the filename.
    """
    name = outline_file.name.lower()
    if "promice" in name:
        return "promice"
    if "simple" in name:
        return "simple"
    if "buffered" in name:
        return "buffered"
    return "detailed"


def default_stride(outline_file: Path) -> int:
    """
    Determine the default stride for the outline file.

    This is needed because the PROMICE outline is very dense, so it will
    inadvertently create very small triangles if not combined into fewer
    segments.

    Parameters
    ----------
    outline_file : Path
        Path to the outline file.

    Returns
    -------
    int
        Default stride for the outline file.
    """
    return {"promice": 1, "simple": 32}.get(outline_kind(outline_file), 2)


def default_simplify_tolerance(outline_file: Path) -> float | None:
    """
    Determine the default simplification tolerance for the outline file.

    Parameters
    ----------
    outline_file : Path
        Path to the outline file.

    Returns
    -------
    float or None
        Default simplification tolerance for the outline file.
    """
    if outline_kind(outline_file) == "promice":
        return PROMICE_SIMPLIFY_TOLERANCE
    return None


def load_velocity(data_dir: Path) -> tuple[xr.DataArray, xr.DataArray]:
    """
    Load the x and y velocity fields from the data directory.

    Parameters
    ----------
    data_dir : Path
        Path to the data directory.

    Returns
    -------
    tuple[xr.DataArray, xr.DataArray]
        x and y velocity fields.
    """
    vx_path = data_dir / MULTIYEAR_VELOCITY_DIR / VX_FILENAME
    vy_path = data_dir / MULTIYEAR_VELOCITY_DIR / VY_FILENAME
    if not vx_path.exists() and (data_dir / VX_FILENAME).exists():
        vx_path = data_dir / VX_FILENAME
    if not vy_path.exists() and (data_dir / VY_FILENAME).exists():
        vy_path = data_dir / VY_FILENAME
    for path in (vx_path, vy_path):
        if not path.exists():
            raise FileNotFoundError(path)

    vx = rxr.open_rasterio(vx_path, masked=True)
    vy = rxr.open_rasterio(vy_path, masked=True)
    return vx, vy


def compute_metric(vx: xr.DataArray, vy: xr.DataArray) -> xr.DataArray:
    """
    Compute the metric for a given velocity field.

    Parameters
    ----------
    vx : xarray.DataArray
        x velocity.
    vy : xarray.DataArray
        y velocity.

    Returns
    -------
    xarray.DataArray
        Metric for the velocity field.
    """
    dudx = vx.differentiate("x")
    dudy = vx.differentiate("y")
    dvdx = vy.differentiate("x")
    dvdy = vy.differentiate("y")

    mag_eps = np.abs(dudx) + np.abs(dvdy) + 0.5 * (np.abs(dudy) + np.abs(dvdx))
    mag_eps = (mag_eps * 20 + np.sqrt(vx**2.0 + vy**2.0) / 2000) ** 0.75
    return mag_eps[0, ::2, ::2].fillna(1.0e-8)


def load_outline(
    outline_file: Path, simplify_tolerance: float | None = None
) -> gpd.GeoDataFrame:
    """
    Load the outline file in the Greenland stereographic projection.

    Parameters
    ----------
    outline_file : Path
        Path to the outline file.
    simplify_tolerance : float or None
        Optional tolerance for simplifying the outline geometry.

    Returns
    -------
    geopandas.GeoDataFrame
        Outline geometry in EPSG:3413.
    """
    if not outline_file.exists():
        raise FileNotFoundError(outline_file)

    outline = gpd.read_file(outline_file).to_crs(EPSG_GREENLAND)
    if simplify_tolerance is not None:
        outline.loc[outline.index[0], "geometry"] = simplify(
            outline.geometry.iloc[0], simplify_tolerance
        )
    return outline


def outline_polygon(outline: gpd.GeoDataFrame) -> Polygon:
    """
    Extract the polygon to mesh from an outline GeoDataFrame.

    Parameters
    ----------
    outline : geopandas.GeoDataFrame
        Outline containing a polygon or multipolygon geometry.

    Returns
    -------
    shapely.geometry.Polygon
        Polygon geometry to mesh. If the outline is a multipolygon, the
        largest polygon part is used.
    """
    geometry = outline.geometry.iloc[0]
    if isinstance(geometry, Polygon):
        return geometry
    if isinstance(geometry, MultiPolygon):
        return max(geometry.geoms, key=lambda polygon: polygon.area)
    raise TypeError(f"Expected Polygon or MultiPolygon, got {geometry.geom_type}")


def densify_closed_line(coords: FloatArray, stride: int) -> FloatArray:
    """
    Densify a closed line by a given stride.

    Needed because one line in the outline may span termini and non-terminus
    segments, so we need to densify it to ensure that we can label these
    differently.

    Parameters
    ----------
    coords : numpy.ndarray
        Coordinates of the line to densify.
    stride : int
        Number of segments to split each edge into.

    Returns
    -------
    numpy.ndarray
        Densified coordinates of the line.
    """
    dense = np.zeros(((coords.shape[0] - 1) * stride + 1, 2))
    dense[0::stride, :] = coords
    for i in range(1, stride):
        dense[i::stride, :] = (
            coords[:-1, :] * (stride - i) + coords[1:, :] * i
        ) / stride
    return dense


def segment_midpoints_and_normals(
    coords: FloatArray,
) -> tuple[xr.DataArray, xr.DataArray, FloatArray, FloatArray]:
    """
    Compute segment midpoints and outward normal vectors.

    Parameters
    ----------
    coords : numpy.ndarray
        Coordinates of the densified outline.

    Returns
    -------
    tuple[xarray.DataArray, xarray.DataArray, numpy.ndarray, numpy.ndarray]
        Midpoint x coordinates, midpoint y coordinates, x normal components,
        and y normal components.
    """
    x = coords[:, 0]
    y = coords[:, 1]
    eval_x = xr.DataArray((x[1:] + x[:-1]) / 2.0, dims="v")
    eval_y = xr.DataArray((y[1:] + y[:-1]) / 2.0, dims="v")

    normal_x = -(y[1:] - y[:-1])
    normal_y = x[1:] - x[:-1]
    norm = np.sqrt(normal_x**2.0 + normal_y**2.0)
    normal_x /= norm
    normal_y /= norm
    return eval_x, eval_y, normal_x, normal_y


def normal_velocity(
    vx: xr.DataArray,
    vy: xr.DataArray,
    eval_x: xr.DataArray,
    eval_y: xr.DataArray,
    normal_x: FloatArray,
    normal_y: FloatArray,
) -> FloatArray:
    """
    Compute the normal velocity at a given set of points.

    Parameters
    ----------
    vx : xarray.DataArray
        Eastward velocity.
    vy : xarray.DataArray
        Northward velocity.
    eval_x : xarray.DataArray
        X coordinates of the points to evaluate the velocity at.
    eval_y : xarray.DataArray
        Y coordinates of the points to evaluate the velocity at.
    normal_x : numpy.ndarray
        X components of the normal vectors at the evaluation points.
    normal_y : numpy.ndarray
        Y components of the normal vectors at the evaluation points.

    Returns
    -------
    numpy.ndarray
        Normal velocities at the evaluation points.
    """
    vx_out = vx.interp(x=eval_x, y=eval_y)
    vy_out = vy.interp(x=eval_x, y=eval_y)
    normal_vel = vx_out * normal_x + vy_out * normal_y
    values = normal_vel.values.flatten()
    values[np.isnan(values)] = 9999
    return values


def split_outlet_boundaries(
    coords: FloatArray, normal_vel: FloatArray, cutoff_vel: float
) -> tuple[list[FloatArray], list[str], list[bool]]:
    """
    Split outline coordinates into outlet and non-outlet boundary segments.

    Parameters
    ----------
    coords : numpy.ndarray
        Densified outline coordinates.
    normal_vel : numpy.ndarray
        Normal velocity evaluated at segment midpoints.
    cutoff_vel : float
        Velocity above which a segment is treated as an outlet.

    Returns
    -------
    tuple[list[numpy.ndarray], list[str], list[bool]]
        Boundary coordinate segments, boundary names, and flags indicating
        whether each segment is an outlet.
    """
    is_outlet = normal_vel[0] > cutoff_vel
    segments = []
    names = []
    outline_lc_map = []
    start_ind = 0

    for i in range(coords.shape[0] - 1):
        is_last_segment = i == coords.shape[0] - 2
        if ((normal_vel[i] > cutoff_vel) != is_outlet) or is_last_segment:
            label = "Outlet" if is_outlet else "Other"
            names.append(f"{label} {len(segments)}")
            segments.append(coords[start_ind:i, :])
            outline_lc_map.append(bool(is_outlet))
            start_ind = i
            is_outlet = normal_vel[i] > cutoff_vel

    return segments, names, outline_lc_map


def write_boundaries(
    boundaries_file: Path, names: Sequence[str], coords: Sequence[FloatArray]
) -> None:
    """
    Write named boundary segments to a geopackage.

    Parameters
    ----------
    boundaries_file : Path
        Output geopackage path.
    names : sequence of str
        Boundary names.
    coords : sequence of numpy.ndarray
        Boundary coordinate segments.
    """
    geometry = [
        LineString(np.vstack((coords[i - 1][-1, :], coords[i])))
        for i in range(len(coords))
    ]
    df = gpd.GeoDataFrame(data={"Name": names}, geometry=geometry, crs=EPSG_GREENLAND)
    if not boundaries_file.exists():
        boundaries_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_file(boundaries_file, driver="GPKG")


def outline_holes(polygon: Polygon) -> list[FloatArray]:
    """
    Extract interior rings from a polygon.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        Polygon geometry.

    Returns
    -------
    list[numpy.ndarray]
        Interior rings with fewer than four points removed.
    """
    cuts = [np.array(geom.xy).T[:-1, :] for geom in polygon.interiors]
    return [cut for cut in cuts if cut.shape[0] > 3]


def make_mesh(
    coords: Sequence[FloatArray],
    names: Sequence[str],
    cuts: Sequence[FloatArray],
    rough_targ: float,
    fine_targ: float,
    outline_lc: Sequence[float],
) -> None:
    """
    Build the active gmsh model from boundary and hole coordinates.

    Parameters
    ----------
    coords : sequence of numpy.ndarray
        Boundary coordinate segments.
    names : sequence of str
        Physical group names for boundary segments.
    cuts : sequence of numpy.ndarray
        Interior hole coordinate rings.
    rough_targ : float
        Coarse characteristic length.
    fine_targ : float
        Fine characteristic length for cuts and outlet boundaries.
    outline_lc : sequence of float
        Characteristic lengths for each boundary segment.
    """
    outline_list = []
    cut_list = []
    point_num = 1
    line_num = 1

    for line_index, line_coords in enumerate(coords):
        first_pt = point_num
        for pt in line_coords:
            gmsh.model.geo.addPoint(pt[0], pt[1], 0, outline_lc[line_index], point_num)
            point_num += 1

        if line_index == 0 and line_index == len(coords) - 1:
            point_tags = np.r_[np.arange(first_pt, point_num), 1]
        elif line_index == 0:
            point_tags = np.arange(first_pt, point_num)
        elif line_index == len(coords) - 1:
            point_tags = np.r_[first_pt - 1, np.arange(first_pt, point_num), 1]
        else:
            point_tags = np.r_[first_pt - 1, np.arange(first_pt, point_num)]

        lines = []
        for i in range(len(point_tags) - 1):
            gmsh.model.geo.addLine(int(point_tags[i]), int(point_tags[i + 1]), line_num)
            lines.append(line_num)
            line_num += 1
        outline_list.append(lines)

    last_outline = line_num - 1

    for cut in cuts:
        first_pt = point_num
        for pt in cut:
            gmsh.model.geo.addPoint(pt[0], pt[1], 0, fine_targ, point_num)
            point_num += 1

        point_tags = np.r_[np.arange(first_pt, point_num), first_pt]
        lines = []
        for i in range(len(point_tags) - 1):
            gmsh.model.geo.addLine(int(point_tags[i]), int(point_tags[i + 1]), line_num)
            lines.append(line_num)
            line_num += 1
        cut_list.append(lines)

    last_cut = line_num - 1
    outline_num = last_cut + 1
    cut_nums = np.arange(outline_num + 1, outline_num + len(cuts) + 1)
    plane_num = last_cut + 1

    gmsh.model.geo.addCurveLoop(np.arange(1, last_outline + 1), outline_num)
    for i, current_cut in enumerate(cut_list):
        gmsh.model.geo.addCurveLoop(current_cut, cut_nums[i])

    gmsh.model.geo.addPlaneSurface([outline_num] + list(cut_nums), plane_num)

    for name, lines in zip(names, outline_list):
        gmsh.model.geo.addPhysicalGroup(1, lines, name=name)
    for i, lines in enumerate(cut_list):
        gmsh.model.geo.addPhysicalGroup(1, lines, name=f"Cut {i}")

    gmsh.model.geo.addPhysicalGroup(2, [plane_num], name="Surf")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)


class Mesh:
    """
    Lightweight view of the current gmsh triangular mesh.

    Attributes
    ----------
    vtags : numpy.ndarray
        gmsh vertex tags.
    vxyz : numpy.ndarray
        Vertex coordinates.
    triangles_tags : numpy.ndarray
        gmsh triangle tags.
    triangles : numpy.ndarray
        Triangle vertex indices into vxyz.
    """

    vtags: IntArray
    vxyz: FloatArray
    triangles_tags: IntArray
    triangles: IntArray

    def __init__(self) -> None:
        """
        Read mesh nodes and triangles from the active gmsh model.
        """
        self.vtags, vxyz, _ = gmsh.model.mesh.getNodes()
        self.vxyz = vxyz.reshape((-1, 3))
        vmap = {tag: i for i, tag in enumerate(self.vtags)}
        self.triangles_tags, evtags = gmsh.model.mesh.getElementsByType(2)
        evid = np.array([vmap[tag] for tag in evtags])
        self.triangles = evid.reshape((self.triangles_tags.shape[-1], -1))


def compute_size_field(
    nodes: FloatArray,
    triangles: IntArray,
    rough_targ: float,
    fine_targ: float,
    meps: xr.DataArray,
) -> FloatArray:
    """
    Compute target element sizes from mesh centroids and the velocity metric.

    Parameters
    ----------
    nodes : numpy.ndarray
        Mesh node coordinates.
    triangles : numpy.ndarray
        Triangle vertex indices into nodes.
    rough_targ : float
        Maximum target element size.
    fine_targ : float
        Minimum target element size.
    meps : xarray.DataArray
        Velocity-derived mesh metric.

    Returns
    -------
    numpy.ndarray
        Target size for each triangle.
    """
    vxyz = nodes[triangles].mean(axis=1)
    x = xr.DataArray(vxyz[:, 0], dims="H")
    y = xr.DataArray(vxyz[:, 1], dims="H")

    tsize = fine_targ / meps.interp(x=x, y=y).values.flatten()
    tsize[tsize < fine_targ] = fine_targ
    tsize[tsize > rough_targ] = rough_targ
    return tsize


def generate_one_mesh(
    fn_base: str,
    coords: Sequence[FloatArray],
    names: Sequence[str],
    cuts: Sequence[FloatArray],
    outline_lc_map: Sequence[bool],
    meps: xr.DataArray,
    rough_targ: float,
    fine_targ: float,
) -> None:
    """
    Generate one two-pass gmsh mesh with a background size field.

    Parameters
    ----------
    fn_base : str
        Output filename stem.
    coords : sequence of numpy.ndarray
        Boundary coordinate segments.
    names : sequence of str
        Physical group names for boundary segments.
    cuts : sequence of numpy.ndarray
        Interior hole coordinate rings.
    outline_lc_map : sequence of bool
        Flags indicating whether each boundary segment uses the fine length.
    meps : xarray.DataArray
        Velocity-derived mesh metric.
    rough_targ : float
        Maximum target element size.
    fine_targ : float
        Minimum target element size.
    """
    gmsh.initialize()
    try:
        outline_lc = [
            fine_targ if is_outlet else rough_targ for is_outlet in outline_lc_map
        ]

        gmsh.model.add(fn_base + "_orig")
        make_mesh(coords, names, cuts, rough_targ, fine_targ, outline_lc)
        gmsh.write(fn_base + "_raw.msh")
        mesh = Mesh()

        sf_ele = compute_size_field(
            mesh.vxyz, mesh.triangles, rough_targ, fine_targ, meps
        )
        sf_view = gmsh.view.add("mesh size field")
        gmsh.view.addModelData(
            sf_view,
            0,
            fn_base + "_orig",
            "ElementData",
            mesh.triangles_tags,
            sf_ele[:, None],
        )
        gmsh.view.write(sf_view, fn_base + "_sf.pos")

        gmsh.model.add(fn_base)
        make_mesh(coords, names, cuts, rough_targ, fine_targ, outline_lc)

        bg_field = gmsh.model.mesh.field.add("PostView")
        gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
        gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(fn_base + ".msh")
    finally:
        gmsh.finalize()

    for suffix in ("_raw.msh", "_sf.pos"):
        path = Path(fn_base + suffix)
        if path.exists():
            path.unlink()


def plot_mesh(fn_base: str, fig_dir: Path) -> None:
    """
    Plot a generated mesh using Firedrake.

    Parameters
    ----------
    fn_base : str
        Mesh filename stem, without the .msh suffix.
    fig_dir : Path
        Directory where the figure will be saved.
    """
    import firedrake
    import matplotlib.pyplot as plt

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig_path = fig_dir / f"{Path(fn_base).name}.png"
    mesh2d = firedrake.Mesh(fn_base + ".msh")
    fig, ax = plt.subplots(figsize=(14, 14))
    firedrake.triplot(mesh2d, axes=ax)
    ax.axis("equal")
    ax.legend()
    fig.savefig(fig_path, dpi=300)


def main() -> None:
    """
    Run mesh generation from command-line arguments.
    """
    args = parse_args()
    data_dir = args.data_dir.expanduser()
    boundaries_file = args.boundaries_file.expanduser()
    outline_file, outline_name = resolve_outline(args, data_dir)
    stride = default_stride(Path(outline_name)) if args.stride is None else args.stride
    output_dir = args.output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    vx, vy = load_velocity(data_dir)
    meps = compute_metric(vx, vy)
    simplify_tolerance = (
        default_simplify_tolerance(Path(outline_name))
        if args.outline_simplify is None
        else args.outline_simplify
    )
    outline = load_outline(outline_file, simplify_tolerance)

    polygon = outline_polygon(outline)
    rough_coords = np.array(polygon.exterior.xy).T[:, :]
    coords = densify_closed_line(rough_coords, stride)
    eval_x, eval_y, normal_x, normal_y = segment_midpoints_and_normals(coords)
    nv = normal_velocity(vx, vy, eval_x, eval_y, normal_x, normal_y)
    segments, names, outline_lc_map = split_outlet_boundaries(
        coords, nv, args.cutoff_velocity
    )
    write_boundaries(boundaries_file, names, segments)
    cuts = outline_holes(polygon)

    fn_template = str(output_dir / f"greenland_{outline_name}_{{:d}}_{{:d}}")
    for lc in [args.min_lc * 2**i for i in range(args.num_levels)]:
        rough_targ, fine_targ = lc * 10, lc
        fn_base = fn_template.format(rough_targ, fine_targ)
        generate_one_mesh(
            fn_base,
            segments,
            names,
            cuts,
            outline_lc_map,
            meps,
            rough_targ,
            fine_targ,
        )

    if not args.no_plot:
        rough_targ, fine_targ = args.plot_level
        plot_mesh(fn_template.format(rough_targ, fine_targ), args.fig_dir.expanduser())


if __name__ == "__main__":
    main()
