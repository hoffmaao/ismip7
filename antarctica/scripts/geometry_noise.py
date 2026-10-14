r"""Cell-scale noise in the thickness change of an unforced forward.

Long DG0 forwards grow flow-aligned cell-scale streaks in dH. This compares
geometry spaces (``ISMIP7_GEOMETRY_SPACE``) on that measure, from MAPs
inverted in each space:

    geometry_noise.py run NAME YEARS DT     # under mpiexec, from the run's env
    geometry_noise.py measure RESULTS NAME  # serial, after the run

``run`` sets the model up from the MAP (``ISMIP7_INVERSION``), writes the
t = 0 state, and steps ``YEARS`` with no forcing, writing a state every year
(``ISMIP7_CHECKPOINT_EVERY_YR=1``, ``ISMIP7_KEEP_CHECKPOINTS=0``) beside the
timeseries CSV, whose columns carry the mass budget and the front outflux.

``measure`` reads those states and reports, per year, over interior ice (cells
holding at least ``--hmin`` m at t = 0 whose three facet neighbours do too):

* ``rms_dh``: the area-weighted rms of the cell-mean ``dH = h(t) - h(0)``;
* ``rms_nbr``: the area-weighted rms of ``dH_K`` minus the mean of its three
  facet neighbours, the neighbour-scale part of ``dH``;
* ``corr_y1``: the correlation of that neighbour-scale part with the same
  quantity after the first year, which is near 1 when the noise is the first
  year's pattern growing in place.

Cell means are the measure in every space, so DG0 and DG1 compare cell for
cell. A DG1 cell's mean is the mean of its three vertex values.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np


def _run(args):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from simulation import RESULTS_DIR, lc, run_simulation, save_model_state, setup_model

    os.environ.setdefault("ISMIP7_CHECKPOINT_EVERY_YR", "1")
    os.environ.setdefault("ISMIP7_KEEP_CHECKPOINTS", "0")
    ctx = setup_model()
    save_model_state(ctx, os.path.join(
        RESULTS_DIR, f"{args.name}_{lc}_t{args.t0:.1f}.h5"), args.t0)
    run_simulation(ctx, experiment_name=args.name, t_start=args.t0,
                   t_end=args.t0 + args.years, dt=args.dt,
                   output_interval=10 ** 9, checkpoint_interval=10 ** 9,
                   forcing_callback=None)


def _cell_state(path):
    r"""Centroids, areas, cell-mean thickness and cell vertex ids of one
    state file, in that file's cell order."""
    import firedrake as fd
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        h = chk.load_function(mesh, "thickness")
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        fd.SpatialCoordinate(mesh)).dat.data_ro.copy()
    area = fd.assemble(fd.TestFunction(Q0) * fd.dx).dat.data_ro.copy()
    # the cell mean: exact for DG0 and for P1 on a triangle
    mean = fd.Function(Q0).project(h).dat.data_ro.copy()
    cells = mesh.coordinates.cell_node_map().values
    dofs = Q0.cell_node_map().values[:, 0]
    verts = np.empty((dofs.size, 3), dtype=np.int64)
    verts[dofs] = cells[:dofs.size]
    return xy, area, mean, verts


def _order(xy):
    return np.lexsort((np.round(xy[:, 1], 3), np.round(xy[:, 0], 3)))


def facet_neighbours(verts):
    r"""For cells given by their vertex ids, an ``(ncell, 3)`` array of the
    cells across each facet, -1 on the boundary."""
    ncell = verts.shape[0]
    edges = np.concatenate([np.sort(verts[:, [i, j]], axis=1)
                            for i, j in ((0, 1), (1, 2), (0, 2))])
    owner = np.tile(np.arange(ncell), 3)
    key = edges[:, 0] * (verts.max() + 1) + edges[:, 1]
    order = np.argsort(key, kind="stable")
    key, owner = key[order], owner[order]
    nbr = -np.ones((ncell, 3), dtype=np.int64)
    fill = np.zeros(ncell, dtype=np.int64)
    pair = np.flatnonzero(key[1:] == key[:-1])
    for a, b in ((owner[pair], owner[pair + 1]), (owner[pair + 1], owner[pair])):
        for cell, other in zip(a, b):
            nbr[cell, fill[cell]] = other
            fill[cell] += 1
    return nbr


def neighbour_residual(values, nbr):
    r"""``values`` minus the mean over each cell's facet neighbours (NaN
    where a cell has fewer than three)."""
    full = np.all(nbr >= 0, axis=1)
    out = np.full(values.shape, np.nan)
    out[full] = values[full] - values[nbr[full]].mean(axis=1)
    return out


def noise_table(states, hmin=1.0):
    r"""Rows ``(year, rms_dh, rms_nbr, corr_y1)`` and the interior cell count.
    ``states`` is ``(years, means, nbr, area)``: the years, the cell-mean
    thickness of each (the first at t = 0), all in one cell order, the facet
    neighbours and the cell areas in that order."""
    years, means, nbr, area = states
    h0 = means[0]
    full = np.all(nbr >= 0, axis=1)
    ice = h0 >= hmin
    inner = full & ice & np.all(ice[np.where(nbr >= 0, nbr, 0)], axis=1)
    w = area[inner]
    rows = []
    first = None
    for year, mean in zip(years[1:], means[1:]):
        dh = mean - h0
        r = neighbour_residual(dh, nbr)[inner]
        if first is None:
            first = r
        corr = float(np.corrcoef(r, first)[0, 1]) if r.std() > 0 else float("nan")
        rows.append((year, float(np.sqrt((w * dh[inner] ** 2).sum() / w.sum())),
                     float(np.sqrt((w * r ** 2).sum() / w.sum())), corr))
    return rows, int(inner.sum())


def _measure(args):
    pat = re.compile(re.escape(f"{args.name}_") + r"\d+_t(-?\d+\.\d+)\.h5$")
    files = sorted(((float(m.group(1)), f) for f in glob.glob(
        os.path.join(args.results, f"{args.name}_*_t*.h5"))
        if (m := pat.search(os.path.basename(f)))))
    if len(files) < 2:
        raise SystemExit(f"need the t=0 state and at least one year of {args.name}")
    years, means = [], []
    nbr = area = ref = None
    for year, path in files:
        xy, a, mean, verts = _cell_state(path)
        order = _order(xy)
        if ref is None:
            ref = xy[order]
            area = a[order]
            inv = np.empty_like(order)
            inv[order] = np.arange(order.size)
            nbr_raw = facet_neighbours(verts)
            nbr = np.where(nbr_raw >= 0, inv[np.maximum(nbr_raw, 0)], -1)[order]
        elif not np.allclose(xy[order], ref, atol=1e-3):
            raise SystemExit(f"{path} is on another mesh")
        years.append(year)
        means.append(mean[order])
    rows, n = noise_table((years, means, nbr, area), hmin=args.hmin)
    print(f"# {args.name}: {n} interior ice cells (h0 >= {args.hmin:g} m)")
    print("year,rms_dh_m,rms_nbr_m,corr_y1")
    for year, rms_dh, rms_nbr, corr in rows:
        print(f"{year:.1f},{rms_dh:.4f},{rms_nbr:.4f},{corr:.4f}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run")
    run.add_argument("name")
    run.add_argument("years", type=float)
    run.add_argument("dt", type=float)
    run.add_argument("--t0", type=float, default=2015.0)
    meas = sub.add_parser("measure")
    meas.add_argument("results")
    meas.add_argument("name")
    meas.add_argument("--hmin", type=float, default=1.0)
    args = parser.parse_args(argv)
    (_run if args.cmd == "run" else _measure)(args)


if __name__ == "__main__":
    main()
