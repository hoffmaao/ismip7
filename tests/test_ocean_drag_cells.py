r"""The floor-cell ocean drag acts only in open water outside the t=0 extent
that no ice cell touches, so no floating ice and no front vertex feels it."""
import numpy as np
import pytest

from icepack2_tools.front import ocean_drag_cells


def _chain_neighbours(n):
    r"""``neighbours_of`` for cells 0..n-1 in a row, each touching the next."""
    def neighbours_of(mask):
        out = np.zeros(n, bool)
        out[1:] |= mask[:-1]
        out[:-1] |= mask[1:]
        return out
    return neighbours_of


def test_no_drag_on_ice_or_the_water_row_the_front_shares():
    ice = np.zeros(10, bool); ice[3:6] = True
    drag = ocean_drag_cells(ice, _chain_neighbours(10), extent0=ice)
    assert list(np.flatnonzero(drag)) == [0, 1, 7, 8, 9]


def test_a_retreat_inside_the_t0_extent_stays_drag_free():
    r"""Cells that held ice at t=0 and have emptied are water now, but the
    front that retreated through them must not meet drag there."""
    extent0 = np.zeros(10, bool); extent0[2:8] = True
    ice = np.zeros(10, bool); ice[2:4] = True
    drag = ocean_drag_cells(ice, _chain_neighbours(10), extent0)
    assert list(np.flatnonzero(drag)) == [0, 8, 9]


def test_ice_that_advances_past_the_t0_extent_carries_no_drag():
    extent0 = np.zeros(10, bool); extent0[2:5] = True
    ice = np.zeros(10, bool); ice[2:7] = True
    drag = ocean_drag_cells(ice, _chain_neighbours(10), extent0)
    assert list(np.flatnonzero(drag)) == [0, 8, 9]


def test_on_a_mesh_the_gate_is_open_water_a_cell_away_from_the_ice():
    r"""Through the parallel-safe facet transfer the run uses: an ice square in
    the middle of a mesh, drag nowhere on it, nowhere on a cell sharing a
    facet with it, and everywhere else."""
    fd = pytest.importorskip("firedrake")
    from icepack2_tools.front import facet_neighbours

    mesh = fd.UnitSquareMesh(12, 12)
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    x = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    ice = (np.abs(x[:, 0] - 0.5) < 0.2) & (np.abs(x[:, 1] - 0.5) < 0.2)
    neighbours_of = facet_neighbours(Q0)
    drag = ocean_drag_cells(ice, neighbours_of, extent0=ice)
    touching = neighbours_of(ice)
    assert ice.any() and (touching & ~ice).any()
    assert not (drag & ice).any()                    # no drag on any ice cell
    assert not (drag & touching).any()               # nor on the row the front shares
    assert np.array_equal(drag, ~ice & ~touching)    # and on every other cell
