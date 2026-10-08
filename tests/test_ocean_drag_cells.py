r"""The floor-cell ocean drag acts only in open water outside the t=0 extent
that no ice cell touches. The default vertex gate keeps it off every node of
the ice; under ISMIP7_DRAG_GATE=facet a water cell touching the ice at one
vertex still carries it, on that front vertex."""
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


def _cell_nodes_touching(mesh, cells):
    import firedrake as fd
    Q1 = fd.FunctionSpace(mesh, "CG", 1)
    cn = Q1.cell_node_map().values
    nodes = np.zeros(Q1.dof_dset.size, bool)
    nodes[np.unique(cn[cells])] = True
    return cn, nodes


def test_the_facet_gate_drags_front_vertices_and_the_vertex_gate_does_not():
    r"""Issue #153: a water cell sharing only a vertex with the ice is outside
    the facet gate, so its drag acts on that front vertex of the CG1
    velocity. On the 20 km buffered 2 km mesh that held RC's shelves 27 %
    below the observed speed. The vertex gate leaves no drag cell holding a
    node of an ice cell, and drags every other open-water cell."""
    fd = pytest.importorskip("firedrake")
    from icepack2_tools.front import facet_neighbours, vertex_neighbours

    mesh = fd.UnitSquareMesh(12, 12)
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    x = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro
    ice = (np.abs(x[:, 0] - 0.5) < 0.2) & (np.abs(x[:, 1] - 0.5) < 0.2)
    cn, ice_nodes = _cell_nodes_touching(mesh, ice)

    facet_drag = ocean_drag_cells(ice, facet_neighbours(Q0), extent0=ice)
    leaking = facet_drag & ice_nodes[cn].any(axis=1)
    assert leaking.any()

    by_vertex = vertex_neighbours(Q0)
    vertex_drag = ocean_drag_cells(ice, by_vertex, extent0=ice)
    assert not (vertex_drag & ice_nodes[cn].any(axis=1)).any()
    assert np.array_equal(vertex_drag, ~ice & ~ice_nodes[cn].any(axis=1))
    assert np.array_equal(by_vertex(ice), ice_nodes[cn].any(axis=1))
    assert (facet_drag & ~vertex_drag).sum() == leaking.sum()


@pytest.mark.parametrize("value, expected",
                         [(None, "vertex"), ("facet", "facet"), (" Vertex ", "vertex")])
def test_the_drag_gate_knob(monkeypatch, value, expected):
    from icepack2_tools.runconfig import drag_gate
    if value is None:
        monkeypatch.delenv("ISMIP7_DRAG_GATE", raising=False)
    else:
        monkeypatch.setenv("ISMIP7_DRAG_GATE", value)
    assert drag_gate() == expected
    monkeypatch.setenv("ISMIP7_DRAG_GATE", "node")
    with pytest.raises(ValueError, match="ISMIP7_DRAG_GATE"):
        drag_gate()


def test_a_forward_runs_the_gate_its_map_records(monkeypatch):
    from icepack2_tools.runconfig import forward_drag_gate
    monkeypatch.delenv("ISMIP7_DRAG_GATE", raising=False)
    assert forward_drag_gate("facet") == "facet"
    assert forward_drag_gate(b"vertex") == "vertex"
    monkeypatch.setenv("ISMIP7_DRAG_GATE", "facet")
    assert forward_drag_gate("facet") == "facet"
    with pytest.raises(RuntimeError, match="follows its MAP"):
        forward_drag_gate("vertex", source="m.h5")
    with pytest.raises(RuntimeError, match="drag_gate='node'"):
        forward_drag_gate("node")


def test_a_map_without_a_gate_runs_the_knob(monkeypatch):
    from icepack2_tools.runconfig import forward_drag_gate
    monkeypatch.delenv("ISMIP7_DRAG_GATE", raising=False)
    # no record, or an inversion that dragged no cell: the default
    assert forward_drag_gate(None) == "vertex"
    assert forward_drag_gate("none") == "vertex"
    # a forward state from before the record ran facet
    assert forward_drag_gate(None, restart=True) == "facet"
    monkeypatch.setenv("ISMIP7_DRAG_GATE", "facet")
    assert forward_drag_gate("none") == "facet"
    monkeypatch.setenv("ISMIP7_DRAG_GATE", "vertex")
    assert forward_drag_gate(None, restart=True) == "vertex"
