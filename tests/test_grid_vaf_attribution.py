r"""grid_vaf_attribution.py: the organisers' volume above flotation on 8 km
pixel means against the model's, on overlap matrices small enough to work by
hand (issue #99).

Each pixel below is covered by two cells of half its area, except the domain
edge pixel, which one cell covers by half. The cases are the three ways the
grid parts from the mesh: a grounding line whose open-water half cancels its
grounded half, a bed that crosses sea level inside a grounded pixel, and a
pixel the mesh only partly covers.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sp = pytest.importorskip("scipy.sparse")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "antarctica" / "scripts"))
import grid_vaf_attribution as gv  # noqa: E402

HALF = gv.PIXEL_AREA / 2


def hf(b):
    return max(-b, 0.0) * 1024.0 / 917.0


def two_halves(cells):
    r"""One pixel per pair of cells, each cell half the pixel."""
    n = len(cells) // 2
    rows = np.repeat(np.arange(n), 2)
    W = sp.csr_matrix((np.full(2 * n, HALF), (rows, np.arange(2 * n))), shape=(n, 2 * n))
    h = np.array([c[0] for c in cells], dtype=float)
    b = np.array([c[1] for c in cells], dtype=float)
    return W, h, b, np.full(2 * n, HALF)


def test_open_water_cancels_the_grounded_half_of_a_grounding_line_pixel():
    W, h, b, area = two_halves([(1000.0, -500.0), (0.0, -1000.0)])
    e = gv.evaluate(W, h, b, area)
    mesh_m = 0.5 * (1000.0 - hf(-500.0))
    # the pixel means, 500 m of ice over a -750 m bed, float: the tool reads none
    assert e["grid_vaf_gt"] == 0.0
    assert e["mesh_vaf_gt"] == pytest.approx(mesh_m * gv.GT_PER_PIXEL_M)
    assert e["gl_pixels"] == 1 and e["edge_pixels"] == 0
    assert e["gl_gt"] == pytest.approx(-mesh_m * gv.GT_PER_PIXEL_M)
    assert e["gl_cancellation_gt"] == pytest.approx(e["gl_gt"])
    assert e["gl_convexity_gt"] == 0.0
    assert e["gl_supply_gt"] == pytest.approx(mesh_m * gv.GT_PER_PIXEL_M)
    assert e["mesh_sum_check"] == pytest.approx(0.0, abs=1e-12)


def test_removing_a_shelf_moves_the_grid_and_leaves_the_mesh():
    # grounded 2000 m on a -500 m bed beside a 600 m shelf over a -600 m bed
    W, h, b, area = two_halves([(2000.0, -500.0), (600.0, -600.0)])
    base = gv.evaluate(W, h, b, area, keep=True)
    cells = base["_cells"]
    assert list(cells["grounded"]) == [True, False]
    states = {name: h2 for name, h2, _ in gv.perturbations(h, b, cells["grounded"], cells["ice"])}
    gone = gv.evaluate(W, states["shelf_removed"], b, area)
    assert gone["mesh_vaf_gt"] == pytest.approx(base["mesh_vaf_gt"])
    # both beds are below sea level, so h_f is linear across the pixel and all
    # of grid minus mesh is cancellation: -35.0 m with the shelf, -335.0 without
    mean_hf = 0.5 * (hf(-500.0) + hf(-600.0))
    mesh_m = 0.5 * (2000.0 - hf(-500.0))
    assert base["gl_gt"] == pytest.approx((1300.0 - mean_hf - mesh_m) * gv.GT_PER_PIXEL_M)
    assert gone["gl_gt"] == pytest.approx((1000.0 - mean_hf - mesh_m) * gv.GT_PER_PIXEL_M)
    assert gone["gl_gt"] - base["gl_gt"] == pytest.approx(-300.0 * gv.GT_PER_PIXEL_M)
    thinned = gv.evaluate(W, states["shelf_thinned_50"], b, area)
    assert thinned["gl_gt"] - base["gl_gt"] == pytest.approx(-150.0 * gv.GT_PER_PIXEL_M)


def test_a_bed_across_sea_level_puts_the_grid_above_the_mesh():
    # fully grounded: 1000 m on a +200 m bed beside 1000 m on a -600 m bed
    W, h, b, area = two_halves([(1000.0, 200.0), (1000.0, -600.0)])
    e = gv.evaluate(W, h, b, area)
    convexity_m = 0.5 * hf(-600.0) - hf(-200.0)
    assert e["gl_pixels"] == 0
    assert e["interior_gt"] == pytest.approx(convexity_m * gv.GT_PER_PIXEL_M)
    assert e["interior_convexity_gt"] == pytest.approx(e["interior_gt"])
    assert e["interior_cancellation_gt"] == pytest.approx(0.0, abs=1e-9)


def test_a_partly_covered_pixel_is_domain_edge():
    # one grounded cell covers half the pixel: lithk is a whole-pixel mean and
    # topg a mean over the covered half, so the diluted thickness floats
    W = sp.csr_matrix(np.array([[HALF]]))
    e = gv.evaluate(W, np.array([1000.0]), np.array([-500.0]), np.array([HALF]))
    assert e["edge_pixels"] == 1 and e["gl_pixels"] == 0
    assert e["edge_gt"] == pytest.approx(-0.5 * (1000.0 - hf(-500.0)) * gv.GT_PER_PIXEL_M)


def test_a_retreat_takes_marine_ice_near_flotation_and_the_shelves():
    # marine 10 m above flotation, marine 300 m above, 20 m of ice on land, a shelf
    h = np.array([hf(-500.0) + 10.0, hf(-500.0) + 300.0, 20.0, 300.0])
    b = np.array([-500.0, -500.0, 100.0, -600.0])
    grounded = np.array([True, True, True, False])
    ice = np.ones(4, dtype=bool)
    states = {name: (h2, gone) for name, h2, gone in gv.perturbations(h, b, grounded, ice)}
    h25, gone25 = states["retreat_haf25"]
    assert list(gone25) == [True, False, False, False]
    assert list(h25) == [0.0, h[1], 20.0, 0.0]
    assert list(states["retreat_haf800"][1]) == [True, True, False, False]
    assert not states["shelf_removed"][1].any()
