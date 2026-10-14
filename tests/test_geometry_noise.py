r"""The neighbour-scale noise measure of ``antarctica/scripts/geometry_noise.py``.

Pure numpy on a structured triangulation: the facet neighbours it finds, a
smooth field carrying no neighbour-scale part, and a cell-scale pattern
growing in place scoring a correlation of 1 with its first year.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "antarctica" / "scripts"))
from geometry_noise import facet_neighbours, neighbour_residual, noise_table  # noqa: E402


def _triangulation(n):
    r"""An n x n grid of unit squares, each cut into two triangles: vertex
    ids per cell and the cell centroids."""
    vid = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    cells, xy = [], []
    for i in range(n):
        for j in range(n):
            a, b, c, d = vid[i, j], vid[i + 1, j], vid[i, j + 1], vid[i + 1, j + 1]
            cells += [(a, b, d), (a, d, c)]
            xy += [(i + 2 / 3, j + 1 / 3), (i + 1 / 3, j + 2 / 3)]
    return np.array(cells), np.array(xy)


def test_each_cell_finds_the_cells_across_its_facets():
    cells, _ = _triangulation(3)
    nbr = facet_neighbours(cells)
    for k, row in enumerate(nbr):
        for other in row[row >= 0]:
            assert len(set(cells[k]) & set(cells[other])) == 2
            assert k in nbr[other]
    # interior cells have three neighbours, and the 12 boundary facets one side
    assert (nbr < 0).sum() == 12
    assert np.all((nbr >= 0).sum(axis=1) >= 1)


def test_a_linear_field_has_no_neighbour_scale_part():
    r"""On this mesh the three facet neighbours of a cell average to its own
    centroid, so a linear field is its neighbours' mean exactly."""
    cells, xy = _triangulation(8)
    nbr = facet_neighbours(cells)
    r = neighbour_residual(3.0 + 0.7 * xy[:, 0] - 0.2 * xy[:, 1], nbr)
    full = ~np.isnan(r)
    assert full.any() and np.allclose(r[full], 0.0)


def test_a_cell_scale_pattern_growing_in_place_correlates_with_its_first_year():
    cells, xy = _triangulation(10)
    nbr = facet_neighbours(cells)
    area = np.full(len(cells), 0.5)
    h0 = np.full(len(cells), 500.0)
    rng = np.random.default_rng(3)
    pattern = rng.normal(size=len(cells))
    smooth = 0.1 * xy[:, 0]
    years = [0.0, 1.0, 2.0, 3.0]
    means = [h0 + t * (pattern + smooth) for t in years]
    rows, n = noise_table((years, means, nbr, area))
    assert n > 0
    rms = [r[2] for r in rows]
    np.testing.assert_allclose(rms, [rms[0], 2 * rms[0], 3 * rms[0]], rtol=1e-12)
    np.testing.assert_allclose([r[3] for r in rows], 1.0, rtol=1e-12)
    # the rms of dH counts the smooth part too, the neighbour-scale part not
    bare, _ = noise_table((years, [h0 + t * pattern for t in years], nbr, area))
    np.testing.assert_allclose([r[2] for r in bare], rms, rtol=1e-9)
    assert rows[0][1] > bare[0][1]
