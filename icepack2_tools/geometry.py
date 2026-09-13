r"""Geometry representation helpers shared by the inversion and the forward.

The prognostic geometry (``h``, ``s``, ``b``) lives in the space selected by
``ISMIP7_GEOMETRY_SPACE`` - DG0 by default, so that the momentum solve and the
mass transport use one thickness field and the calving-terminus traction cannot
disagree with the boundary flux. See ``GEOMETRY_DISCRETIZATION.md``.

Three operations need care under DG0 and are collected here so the inversion
and the forward cannot drift apart:

* :func:`sample_to_geometry` - get a raster onto the geometry space as a CELL
  AVERAGE rather than a centroid point sample.
* :func:`cg1_lift` - a bounded, volume-preserving CG1 reconstruction, for the
  few places that genuinely need a pointwise gradient of a cell-wise field.
* :func:`surface_slope` - ``grad(s)`` that works for CG1 or DG0 surfaces.
"""

import weakref

import numpy as np
from mpi4py import MPI

from firedrake import (
    Constant,
    Function,
    FunctionSpace,
    TestFunction,
    assemble,
    dx,
    grad,
    max_value,
)


_LUMPED_MASS = weakref.WeakKeyDictionary()


def _lumped_mass(mesh, Q_cg):
    r"""Lumped CG1 mass vector for this mesh, assembled once.

    It depends only on the mesh, and cg1_lift now runs once per forcing
    callback (~1 per timestep), so re-assembling it every call is half the
    per-call cost for nothing. Keyed on the MESH, which lives for the whole
    run: the CG1 space object cg1_lift builds is transient (nothing outside
    holds it), so keying on that would evict the entry on every call. The
    value is the plain array, not the assembled Cofunction, which would keep
    a strong reference back to its own key and pin the mesh forever.
    """
    cached = _LUMPED_MASS.get(mesh)
    if cached is None:
        cached = assemble(TestFunction(Q_cg) * dx).dat.data_ro.copy()
        _LUMPED_MASS[mesh] = cached
    return cached


def cg1_lift(f):
    r"""Lumped-mass CG1 reconstruction of a cell-wise (DG0) field.

    Each CG1 node takes the area-weighted mean of the adjacent cell values.
    Volume-preserving (``int f dx`` is exact) and a convex combination, so it
    cannot overshoot - unlike an L2 projection, which does (projecting the DG0
    thermomechanical fluidity prior to CG1 produced a NEGATIVE fluidity, min
    -9.88 against a DG0 range of [1.0, 446.7]).

    It is NOT unbiased on the domain boundary, where the stencil is one-sided
    and pulls boundary nodes toward interior values. On an unbuffered mesh the
    boundary is the calving front, and the terminus traction goes as ``h^2``,
    so using this on the thickness inflated the front (105 -> 200 m at 32 km)
    and multiplied the outflux by ~4.7. Keep it out of the momentum residual
    and out of every flux; it is for diagnostics, fixed reference scalings, and
    parameterizations only.
    """
    mesh = f.function_space().mesh()
    Q_cg = FunctionSpace(mesh, "CG", 1)
    lumped = _lumped_mass(mesh, Q_cg)
    rhs = assemble(TestFunction(Q_cg) * f * dx)
    out = Function(Q_cg)
    out.dat.data[:] = rhs.dat.data_ro / lumped
    return out


def surface_slope(s):
    r"""``grad(s)`` valid for a CG1 *or* DG0 surface.

    A DG0 surface has an identically zero cell gradient - UFL folds it away -
    because its slope lives entirely in the inter-cell jumps, which the
    momentum balance picks up weakly through its ``jump(s, nu)`` facet term.
    Callers that need a POINTWISE slope (a friction anchor, a melt
    parameterization) get one from a CG1 reconstruction instead.
    """
    if s.function_space().ufl_element().degree() == 0:
        return grad(cg1_lift(s))
    return grad(s)


def _lattice_bary(n):
    r"""Barycentric centroids of the ``n**2`` equal-area sub-triangles of a
    uniform n-fold split of a triangle, shape ``(n*n, 3)``."""
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    up = (i + j) <= n - 1                     # n(n+1)/2 upright sub-triangles
    dn = (i + j) <= n - 2                     # n(n-1)/2 inverted ones
    a = np.concatenate([(i[up] + 1.0 / 3.0) / n, (i[dn] + 2.0 / 3.0) / n])
    b = np.concatenate([(j[up] + 1.0 / 3.0) / n, (j[dn] + 2.0 / 3.0) / n])
    return np.column_stack([a, b, 1.0 - a - b])


def raster_cell_mean(dataset, Q_dg, floor=None, nmax=64, chunk=2048):
    r"""Mean of a rasterio raster over each cell of a DG0 space.

    Each owned cell is split into ``n**2`` equal-area sub-triangles with
    ``n = ceil(sqrt(cell_area / pixel_area))`` (capped at ``nmax``), the raster
    is looked up at every sub-triangle centroid by *pixel containment*, and the
    cell value is the mean. Sampling density therefore tracks pixel density: a
    2 km cell on 500 m BedMachine gets 9 samples for its ~7 pixels, a 20 km
    interior cell gets 729 for its ~1600. Contrast :func:`sample_to_geometry`
    with ``method="vertex"``, which reads three pixels per cell whatever its
    size, so the interior's sub-grid roughness aliases straight into the DG0
    facet jumps that ARE the driving stress.

    One raster window covering the rank's cells is read. With a locality-
    preserving partition that window is small; with PETSc's ``simple``
    partitioner every rank's cells are scattered continent-wide and the window
    is the whole raster (711 MB float32 for BedMachine), which is tolerated
    rather than fixed here.

    Masked / nodata pixels are excluded from the mean; a cell with no valid
    sample is left NaN for the caller to fill.
    """
    from rasterio.transform import rowcol
    from rasterio.windows import from_bounds

    mesh = Q_dg.mesh()
    X = mesh.coordinates.dat.data_ro_with_halos[:, :2]
    cells = mesh.coordinates.cell_node_map().values          # owned cells -> nodes
    dofs = Q_dg.cell_node_map().values[:, 0]                 # owned cells -> dof
    out = Function(Q_dg)
    nc = cells.shape[0]
    if nc == 0:
        return out
    tri = X[cells]                                           # (nc, 3, 2)
    v0, v1, v2 = tri[:, 0], tri[:, 1], tri[:, 2]
    area = 0.5 * np.abs((v1[:, 0] - v0[:, 0]) * (v2[:, 1] - v0[:, 1])
                        - (v2[:, 0] - v0[:, 0]) * (v1[:, 1] - v0[:, 1]))
    apix = abs(dataset.res[0] * dataset.res[1])
    n = np.clip(np.ceil(np.sqrt(area / apix)).astype(int), 1, nmax)

    flat = tri.reshape(-1, 2)
    bnd = dataset.bounds
    rx, ry = abs(dataset.res[0]), abs(dataset.res[1])
    win = from_bounds(max(flat[:, 0].min() - 2 * rx, bnd.left),
                      max(flat[:, 1].min() - 2 * ry, bnd.bottom),
                      min(flat[:, 0].max() + 2 * rx, bnd.right),
                      min(flat[:, 1].max() + 2 * ry, bnd.top),
                      transform=dataset.transform)
    win = win.round_lengths(op="ceil").round_offsets(op="floor")
    arr = dataset.read(1, window=win, masked=True)
    arr = np.ma.filled(arr.astype("f4"), np.nan)
    row0, col0 = int(win.row_off), int(win.col_off)
    H, W = arr.shape

    means = np.full(nc, np.nan)
    for nn in np.unique(n):
        sel = np.flatnonzero(n == nn)
        lam = _lattice_bary(int(nn))                         # (k, 3)
        for s0 in range(0, sel.size, chunk):
            idx = sel[s0:s0 + chunk]
            P = np.einsum("kb,cbd->ckd", lam, tri[idx])      # (c, k, 2)
            rows, cols = rowcol(dataset.transform,
                                P[..., 0].ravel(), P[..., 1].ravel())
            r = np.clip(np.asarray(rows) - row0, 0, H - 1)
            c = np.clip(np.asarray(cols) - col0, 0, W - 1)
            vals = arr[r, c].reshape(idx.size, -1)
            with np.errstate(all="ignore"):
                means[idx] = np.nanmean(vals, axis=1, dtype="f8")
    if floor is not None:
        means = np.maximum(means, floor)
    out.dat.data[dofs] = means
    return out


def sample_to_geometry(raster, Q_g, Q_cg, floor=None, method="vertex"):
    r"""Sample a raster onto the geometry space ``Q_g`` as a cell average.

    ``raster`` is either a rasterio dataset, or (legacy) a closure
    ``raster_fn(space)`` returning the raster interpolated onto ``space``.

    ``method`` selects the DG0 cell value (``runconfig.RASTER_SAMPLES``):
    ``"vertex"`` projects the CG1 vertex interpolant (three pixels per cell);
    ``"cell_mean"`` is :func:`raster_cell_mean`, the raster's mean over the
    cell. Under CG1 geometry there is no cell, so ``method`` is ignored.

    For CG1 geometry this is just the nodal interpolant. For DG0 it is NOT the
    obvious ``icepack.interpolate(raster, Q_g)``: a DG0 dof sits at the cell
    centroid, so that would take a ONE-POINT sample of a 500 m BedMachine
    raster per (at 32 km) 32 km cell. Since the DG0 driving stress is entirely
    the facet jump in ``s``, that sampling noise is read as slope. Measured at
    32 km, centroid sampling against the cell average:

        rms |jump s|     366 m  vs  257 m       (42% rougher)
        peakedness       2.04   vs  1.49
        front <h>        209 m  vs  153 m       (BedMachine ice front: 145
                                                 all / 167 floating, median 152)
        |driving force|  7.7%   vs  1.0% from the CG1 value

    and the rough version failed to converge in 200 Newton iterations. The L2
    projection of the CG1 interpolant IS the cell average, which is what a DG0
    field means, so use that.

    ``floor`` optionally clamps the field from below (thickness >= h_clamp).
    The two paths apply it at opposite ends of the averaging and so disagree
    when ``floor > 0``: ``"vertex"`` clamps the CG1 interpolant and then
    projects (clamp, then average), while ``"cell_mean"`` averages the raster
    over the cell and clamps that average (average, then clamp). Only the
    second can return exactly ``floor``; the first returns it only where every
    vertex is at or below the floor. They coincide wherever the raster already
    exceeds the floor, which is why the inversion sees no difference
    (``ISMIP7_H_CLAMP`` defaults to 0 and BedMachine thickness is
    non-negative); the gap is reachable from the budd_legacy cold start, whose
    ``h_clamp_init`` is 10 m.
    """
    if callable(raster):
        raster_fn, dataset = raster, None
    else:
        import icepack
        dataset = raster
        raster_fn = lambda sp: icepack.interpolate(dataset, sp)  # noqa: E731

    is_dg0 = Q_g.ufl_element() != Q_cg.ufl_element()
    if is_dg0 and method == "cell_mean":
        if dataset is None:
            raise TypeError("method='cell_mean' needs a rasterio dataset, "
                            "not a closure")
        out = raster_cell_mean(dataset, Q_g, floor=floor)
        bad = np.isnan(out.dat.data_ro)
        # The fill is a collective L2 projection, so the branch must be taken
        # by every rank or none: a rank-local `bad.any()` hangs under MPI when
        # the nodata cells (or the owned cells) are not spread over all ranks.
        if Q_g.mesh().comm.allreduce(bool(bad.any()), op=MPI.LOR):
            # No valid pixel under the cell (nodata): take the vertex value.
            fill = Function(Q_g).project(raster_fn(Q_cg))
            out.dat.data[bad] = fill.dat.data_ro[bad]
            if floor is not None:
                out.dat.data[bad] = np.maximum(out.dat.data_ro[bad], floor)
        return out
    if method not in ("vertex", "cell_mean"):
        raise ValueError(f"unknown raster sampling method {method!r}")

    field = raster_fn(Q_cg)
    if floor is not None:
        field = Function(Q_cg).interpolate(max_value(field, Constant(floor)))
    if not is_dg0:
        return field if isinstance(field, Function) else Function(Q_cg).interpolate(field)
    return Function(Q_g).project(field)
