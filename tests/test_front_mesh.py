"""A buffered mesh that follows a marine ice front (issue #167).

The front is taken from a classified ice mask (obs_icemask.classify), cut into
the runs whose ice-free side is marine, and embedded in the gmsh surface, so
the mesh's nodes and edges lie on it and no cell holds ice on one side of the
front and ocean on the other.
"""

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "antarctica", "scripts"))

gmsh = pytest.importorskip("gmsh")
pytest.importorskip("rasterio")
pytest.importorskip("shapely")

from rasterio.transform import Affine  # noqa: E402
from shapely.geometry import LinearRing, Point, Polygon, box  # noqa: E402

from icepack2_tools import mesh as M  # noqa: E402
from icepack2_tools.naming import (  # noqa: E402
    mesh_basename, mesh_front, parse_mesh_basename,
)
from icepack2_tools.obs_icemask import ICE, LAND, MARINE, classify  # noqa: E402

PX = 500.0
N = 200                                        # a 100 km square of pixels
TRANSFORM = Affine(PX, 0.0, -N * PX / 2, 0.0, -PX, N * PX / 2)
OUTLINE = box(-48e3, -48e3, 48e3, 48e3)
RADIUS = 30e3


def _centres():
    c = (np.arange(N) + 0.5) * PX - N * PX / 2
    return np.meshgrid(c, -c)                 # rows run north to south


def _disk_classes(land_west):
    x, y = _centres()
    out = np.where(x > 0, MARINE, LAND if land_west else MARINE).astype("i1")
    out[np.hypot(x, y) <= RADIUS] = ICE
    return out


def test_classify_splits_no_ice_into_marine_and_land():
    ice = np.array([1, 0, 0, 0, 0, 0], bool)
    bm = np.array([2, 0, 1, 3, 2, 2])
    bed = np.array([100.0, -500.0, 50.0, -300.0, -20.0, 20.0])
    # ice | ocean | ice-free land | BedMachine shelf | grounded below sea
    # level | grounded above it
    assert classify(ice, bm, bed).tolist() == [ICE, MARINE, LAND, MARINE, MARINE, LAND]


def test_smooth_marine_bridges_short_gaps_and_drops_short_runs():
    flags = np.array([1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0], bool)
    length = np.ones(len(flags))
    out = M.smooth_marine(flags, length, gap=2.0, min_run=3.0)
    # the one-edge land gap is bridged, the long land run stays, and the
    # one-edge marine run inside it goes
    assert out.tolist() == [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]


def test_a_front_half_on_land_is_one_open_run_on_the_marine_side():
    curves, stats = M.extract_marine_front(
        _disk_classes(land_west=True), TRANSFORM, OUTLINE, 2000.0)
    assert len(curves) == 1 and not curves[0]["closed"]
    xy = curves[0]["xy"]
    assert xy[:, 0].min() > -2 * 2000.0
    r = np.hypot(xy[:, 0], xy[:, 1])
    assert np.all(np.abs(r - RADIUS) < 2 * PX)
    assert stats["length_km"] == pytest.approx(np.pi * RADIUS / 1e3, rel=0.06)
    gaps = np.hypot(*np.diff(xy, axis=0).T)
    assert gaps.max() <= 2000.0 + 1e-6


def test_small_islands_and_holes_are_cleaned_away():
    classes = _disk_classes(land_west=False)
    x, y = _centres()
    classes[np.hypot(x - 40e3, y + 40e3) <= 1.5e3] = ICE       # 7 km^2 island
    classes[np.hypot(x, y) <= 1.5e3] = MARINE                  # 7 km^2 hole
    curves, _ = M.extract_marine_front(classes, TRANSFORM, OUTLINE, 2000.0)
    assert len(curves) == 1 and curves[0]["closed"]


def _mesh_with_front(curves, size):
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    try:
        gmsh.model.add("front")
        corners = list(OUTLINE.exterior.coords)[:-1]
        for i, (x, y) in enumerate(corners, start=1):
            gmsh.model.geo.addPoint(x, y, 0.0, size, i)
        for i in range(1, 5):
            gmsh.model.geo.addLine(i, i % 4 + 1, i)
        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 5)
        gmsh.model.geo.addPlaneSurface([5], 6)
        gmsh.model.geo.synchronize()
        M.embed_front_in_gmsh(curves, 6, 5, 6, size)
        gmsh.option.setNumber("Mesh.MeshSizeMax", size)
        gmsh.model.mesh.generate(2)
        _, xyz, _ = gmsh.model.mesh.getNodes()
        tags, tri = gmsh.model.mesh.getElementsByType(2)
        node_tags, _, _ = gmsh.model.mesh.getNodes()
    finally:
        gmsh.finalize()
    xyz = np.asarray(xyz).reshape(-1, 3)[:, :2]
    index = {int(t): i for i, t in enumerate(node_tags)}
    tri = np.array([index[int(t)] for t in tri]).reshape(-1, 3)
    return xyz, tri


def test_the_mesh_puts_nodes_on_the_front_and_no_cell_straddles_it():
    curves, _ = M.extract_marine_front(
        _disk_classes(land_west=False), TRANSFORM, OUTLINE, 2000.0)
    assert len(curves) == 1 and curves[0]["closed"]
    xyz, tri = _mesh_with_front(curves, 2000.0)
    assert M.front_points_off_mesh(curves, xyz) == 0
    ring = LinearRing(curves[0]["xy"])
    ice = Polygon(ring)
    on = np.array([ring.distance(Point(p)) < 1e-6 for p in xyz])
    inside = np.array([ice.contains(Point(p)) for p in xyz]) & ~on
    outside = ~inside & ~on
    straddles = inside[tri].any(axis=1) & outside[tri].any(axis=1)
    assert not straddles.any()
    assert inside.any() and outside.any()


def test_save_and_load_front_round_trip(tmp_path):
    curves, stats = M.extract_marine_front(
        _disk_classes(land_west=True), TRANSFORM, OUTLINE, 2000.0)
    path = tmp_path / "front.npz"
    M.save_front(curves, path, stats)
    back = M.load_front(path)
    assert len(back) == len(curves)
    for a, b in zip(curves, back):
        assert a["closed"] == b["closed"]
        np.testing.assert_array_equal(a["xy"], b["xy"])


def test_front_meshes_are_named_by_their_edge(monkeypatch):
    name = mesh_basename(10000, 1000, 20000.0, "bm")
    assert name == "antarctica_10000_1000_buffered20000_frontbm"
    assert parse_mesh_basename(name + ".msh") == (10000, 1000, 20000)
    assert mesh_front(name) == "bm"
    assert mesh_front("antarctica_10000_1000_buffered20000") is None
    assert mesh_basename(10000, 1000, 20000.0) == "antarctica_10000_1000_buffered20000"

    from icepack2_tools.runconfig import mesh_front as knob
    from mesh_naming import bndids_filename, mesh_filename
    monkeypatch.delenv("ISMIP7_MESH_FRONT", raising=False)
    assert knob() is None
    assert mesh_filename(10000, 1000, 20000.0).endswith("buffered20000.msh")
    monkeypatch.setenv("ISMIP7_MESH_FRONT", "bm")
    assert knob() == "bm"
    assert mesh_filename(10000, 1000, 20000.0).endswith("buffered20000_frontbm.msh")
    assert bndids_filename(10000, 1000, 20000.0).endswith(
        "boundary_ids_antarctica_10000_1000_buffered20000_frontbm.json")
    assert mesh_filename(10000, 1000, 20000.0, front=None).endswith("buffered20000.msh")
    monkeypatch.setenv("ISMIP7_MESH_FRONT", "2015")
    with pytest.raises(ValueError):
        knob()
