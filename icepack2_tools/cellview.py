r"""Per-cell access to the geometry, for the rules that decide cell by cell.

The front, the forcing mask, the drag gate and the removals all decide per
cell and book ``value x cell area``. Under DG0 a cell's thickness is its one
dof. Under DG1 (``ISMIP7_GEOMETRY_SPACE=dg1``) a cell holds a linear field on
three vertex nodes, and these rules read its CELL MEAN, the mean of its three
vertex values, and act on all three dofs at once. Under DG0 every operation
here is the array arithmetic the DG0 code always did, on the same arrays, so
the DG0 path is unchanged.

The positivity operator :func:`clamp_cells` is the only limiter the DG1
transport applies: a cell whose mean falls below the floor is set to the
floor, booked like the DG0 floor, and a cell whose mean is above it but whose
lowest vertex is below has its slope scaled down until that vertex sits on
the floor (Zhang and Shu 2010), which keeps the mean. A cell already at or
above the floor everywhere is returned unchanged, so the operator never acts
on the thickness extrema a bed trough puts under a smooth surface, which a
vertex-based limiter clips (GEOMETRY_DISCRETIZATION.md, "DG1 geometry").
"""

import numpy as np
from firedrake import Function, FunctionSpace

__all__ = ["CellView", "clamp_cells", "geometry_degree"]


def geometry_degree(space):
    r"""The polynomial degree of a geometry function space: 0 or 1. A degree
    1 space must be the vertex-node DG1 of
    :func:`icepack2_tools.geometry.geometry_function_space`."""
    element = space.ufl_element()
    degree = element.degree()
    if degree == 1 and element.family() == "Discontinuous Lagrange" \
            and element.variant() != "equispaced":
        raise ValueError(
            "DG1 geometry needs its nodes on the vertices "
            "(variant='equispaced'): Firedrake's default DG1 nodes sit inside "
            "the cell, where a nodal floor does not bound the vertex values")
    return degree


class CellView:
    r"""The cells of a geometry space ``Q_g``, indexed like a DG0 array.

    Every per-cell array of the forward is indexed by the owned dofs of a DG0
    space on the mesh (``cell_area``, the front masks, the bed). ``means``
    returns a field in that indexing, and ``set_cells``, ``add_cells`` and
    ``write_nodes`` write through it. Under DG0 the view is the identity on
    ``f.dat.data``. Rank-local: owned cells only, no collectives.
    """

    def __init__(self, Q_g, Q0=None):
        if Q_g.ufl_element().family() != "Discontinuous Lagrange":
            raise ValueError("a cell view needs cell-wise (DG0 or DG1) "
                             "geometry; this space is continuous")
        self.Q_g = Q_g
        self.degree = geometry_degree(Q_g)
        mesh = Q_g.mesh()
        if self.degree == 0:
            self.Q0 = Q_g if Q0 is None else Q0
            self.dofs = None
            return
        self.Q0 = FunctionSpace(mesh, "DG", 0) if Q0 is None else Q0
        n0 = self.Q0.dof_dset.size
        cell_dof0 = self.Q0.cell_node_map().values[:, 0]
        cell_dofs = Q_g.cell_node_map().values
        dofs = np.empty((n0, cell_dofs.shape[1]), dtype=cell_dofs.dtype)
        dofs[cell_dof0] = cell_dofs[:cell_dof0.size]
        self.dofs = dofs

    def nodes(self, f):
        r"""The owned values of ``f`` as ``(ncell, nnode)``, a copy."""
        if self.dofs is None:
            return f.dat.data_ro.copy()[:, None]
        return f.dat.data_ro[self.dofs]

    def means(self, f):
        r"""The cell means of ``f``, a read-only array in DG0 indexing. For
        vertex-node P1 on a triangle the mean is the mean of the vertex
        values. Under DG0 it is ``f.dat.data_ro`` itself."""
        if self.dofs is None:
            return f.dat.data_ro
        return f.dat.data_ro[self.dofs].mean(axis=1)

    def mean_function(self, f, out=None):
        r"""The cell means of ``f`` as a DG0 Function (``out`` if given).
        Under DG0, ``f`` itself."""
        if self.dofs is None:
            return f
        out = Function(self.Q0) if out is None else out
        out.dat.data[:] = self.means(f)
        return out

    def set_cells(self, f, mask, values):
        r"""Set every dof of the cells in ``mask`` to the cell's entry of
        ``values`` (an array in DG0 indexing, or a scalar)."""
        values = np.broadcast_to(np.asarray(values, dtype=float), mask.shape)
        if self.dofs is None:
            f.dat.data[mask] = values[mask]
            return
        f.dat.data[self.dofs[mask]] = values[mask][:, None]

    def add_cells(self, f, delta):
        r"""Add each cell's ``delta`` to every dof of the cell, which shifts
        its mean by ``delta`` and leaves its slope alone."""
        if self.dofs is None:
            f.dat.data[:] += delta
            return
        data = f.dat.data
        data[self.dofs] += np.asarray(delta, dtype=float)[:, None]

    def write_nodes(self, f, values):
        r"""Write ``(ncell, nnode)`` values back through the view."""
        if self.dofs is None:
            f.dat.data[:] = values[:, 0]
            return
        f.dat.data[self.dofs] = values


def clamp_cells(v, floor, *ice_free):
    r"""Floor the cell-wise linear thickness ``v``, shape ``(ncell, nnode)``
    (its vertex values), to ``floor`` in place, except where there is no ice.

    Each ``ice_free`` argument is a boolean mask (or ``None``) of cells that
    hold no ice, floored at 0 instead, exactly as
    :func:`icepack2_tools.front.clamp_thickness` floors a DG0 cell. Per cell:

    * a mean below its floor sets the whole cell to the floor (mass added,
      which the caller measures by integrating before and after);
    * otherwise a vertex below the floor scales the slope by
      ``(mean - floor) / (mean - lowest vertex)``, keeping the mean;
    * a cell at or above the floor everywhere is left as it is.

    With one node per cell this is ``clamp_thickness``. Returns ``v``.
    """
    floors = np.full(v.shape[0], float(floor))
    for mask in ice_free:
        if mask is not None:
            floors[np.asarray(mask, dtype=bool)] = 0.0
    mean = v.mean(axis=1)
    low_mean = mean < floors
    v[low_mean] = floors[low_mean][:, None]
    lowest = v.min(axis=1)
    steep = (lowest < floors) & ~low_mean
    if steep.any():
        m, f = mean[steep], floors[steep]
        theta = (m - f) / (m - lowest[steep])
        scaled = m[:, None] + theta[:, None] * (v[steep] - m[:, None])
        # roundoff can leave a vertex a few ulp under the floor
        v[steep] = np.maximum(scaled, f[:, None])
    return v
