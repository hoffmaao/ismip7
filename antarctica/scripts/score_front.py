#!/usr/bin/env python3
r"""Score a forward run's calving front against the Greene et al. (2022) record.

Read-only. For every state checkpoint named, the run's level set (a signed
distance to the ice front, negative in ice) is compared with the signed
distance to the Greene mask of the epoch nearest the checkpoint's year (the
masks are in the ISMIP7 observations kit, 24 epochs 1997-2021 on a 500 m
grid). The residual

    r = phi_obs - phi_model

is read on the model's front cells, weighted by the front length each cell
carries: r > 0 where the model front lies seaward of the observed one (the
model keeps ice the record says is gone), r < 0 where it lies inland. Its
rms is the distance between the two fronts, its bias their net offset, and
the seaward and inland shares say whether a small rms hides a front that is
right on average and wrong everywhere. The same statistics are given over
the band within ``--band`` of the observed front, weighted by cell area,
and per Mouginot basin. A member of a calving campaign is judged on the mean
front rms over the epochs its run covers, against the pinned-front and
no-calving baselines.

    mpiexec -n 8 python antarctica/scripts/score_front.py --out score.json \
        results/hist_cesm2_waccm_<tag>_2000_t2008.9.h5 ... --plot
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from icepack2_tools.runconfig import obs_data_root  # noqa: E402

DEFAULT_KIT = os.path.join(obs_data_root(), "..", "..", "ISMIP7", "AIS", "obs", "mipkit",
                           "AntarcticaObsISMIP7-v1.2.nc")
DEFAULT_BAND = 10e3


def year_of_days(days):
    import datetime
    t0 = datetime.date(1900, 1, 1)
    out = []
    for d in np.asarray(days):
        day = t0 + datetime.timedelta(days=int(d))
        start = datetime.date(day.year, 1, 1)
        length = (datetime.date(day.year + 1, 1, 1) - start).days
        out.append(day.year + (day - start).days / length)
    return np.array(out)


def stats(r, w):
    total = float(w.sum())
    if not total > 0:
        return {"n": 0}
    seaward, inland = r > 0, r < 0
    return {
        "n": int(len(r)),
        "rms_m": float(np.sqrt((r ** 2 * w).sum() / total)),
        "bias_m": float((r * w).sum() / total),
        "mae_m": float((np.abs(r) * w).sum() / total),
        "seaward_share": float(w[seaward].sum() / total),
        "inland_share": float(w[inland].sum() / total),
        "within_1km": float(w[np.abs(r) < 1e3].sum() / total),
        "within_2km": float(w[np.abs(r) < 2e3].sum() / total),
    }


class Record:
    r"""The Greene masks of the kit, read on rank 0, as signed distances on
    demand (negative in ice), and the Mouginot basins."""

    def __init__(self, path, comm):
        self.comm = comm
        self.kit = None
        if comm.rank == 0:
            import netCDF4
            self.kit = netCDF4.Dataset(path)
            self.x = np.asarray(self.kit["x"][:], float)
            self.y = np.asarray(self.kit["y"][:], float)
            self.years = year_of_days(self.kit["greene_mask_time"][:])
            self.basins = np.asarray(self.kit["mouginot_basins"][:])
        self.years = comm.bcast(self.years if comm.rank == 0 else None, root=0)
        self._cache = {}

    def epoch_near(self, year, tolerance):
        k = int(np.argmin(np.abs(self.years - year)))
        return k if abs(self.years[k] - year) <= tolerance else None

    def signed_distance(self, k):
        r"""Rank 0 only: the exact Euclidean signed distance [m] to the edge
        of the mask of epoch ``k`` on the kit grid."""
        if k not in self._cache:
            from scipy.ndimage import distance_transform_edt
            ice = np.asarray(self.kit["icemask_greene"][k, :, :]) > 0
            dx = float(self.x[1] - self.x[0])
            self._cache[k] = np.where(ice, -distance_transform_edt(ice) * dx,
                                      distance_transform_edt(~ice) * dx).astype(np.float32)
        return self._cache[k]

    def sample(self, field, xy):
        r"""Rank 0 only: nearest-grid-point sample of a kit field at points [m]."""
        ix = np.clip(np.rint((xy[:, 0] - self.x[0]) / (self.x[1] - self.x[0])).astype(int), 0, len(self.x) - 1)
        iy = np.clip(np.rint((xy[:, 1] - self.y[0]) / (self.y[1] - self.y[0])).astype(int), 0, len(self.y) - 1)
        return field[iy, ix]


def model_front(mesh, phi):
    r"""Per-cell front length [m] of the level set: each interior facet
    across which ``phi`` changes sign, counted for both cells."""
    from firedrake import TestFunction, assemble, conditional, dS, lt
    Q0 = phi.function_space()
    v = TestFunction(Q0)
    crossing = conditional(lt(phi("+") * phi("-"), 0.0), 1.0, 0.0)
    # the facet integral of the DG0 test function is the facet's length
    form = crossing * (v("+") + v("-")) * dS
    return assemble(form).dat.data_ro.copy()


def score_state(path, record, band, tolerance, comm):
    import firedrake as fd
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        t = float(chk.get_attr("/", "t_yr"))
        phi = chk.load_function(mesh, name="levelset")
    k = record.epoch_near(t, tolerance)
    if k is None:
        return {"file": os.path.basename(path), "t_yr": t, "epoch": None} if comm.rank == 0 else None
    Q0 = phi.function_space()
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro
    area = fd.assemble(fd.TestFunction(Q0) * fd.dx).dat.data_ro
    front_len = model_front(mesh, phi)
    # every rank's cells to rank 0, where the record is: the statistics are
    # weighted sums, so the gathered arrays are what a serial run would hold
    parts = comm.gather((phi.dat.data_ro.copy(), front_len, area.copy(), np.asarray(xy).copy()), root=0)
    if comm.rank != 0:
        return None
    pm, fl, ar, pts = (np.concatenate([p[i] for p in parts]) for i in range(4))
    po = record.sample(record.signed_distance(k), pts)
    bs = record.sample(record.basins, pts)
    r = po - pm
    ok = np.isfinite(po)
    front = ok & (fl > 0)
    near = ok & (np.abs(po) < band)
    out = {"file": os.path.basename(path), "t_yr": t, "epoch": int(k), "epoch_year": float(record.years[k]),
           "front_km": float(fl[front].sum() / 1e3),
           "front": stats(r[front], fl[front]), "band": stats(r[near], ar[near]), "basins": {}}
    for b in np.unique(bs[front]):
        m = front & (bs == b)
        if m.sum() >= 5:
            out["basins"][str(int(b))] = dict(stats(r[m], fl[m]), front_km=float(fl[m].sum() / 1e3))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("states", nargs="+", help="state checkpoints of one run (any years)")
    ap.add_argument("--kit", default=DEFAULT_KIT)
    ap.add_argument("--band", type=float, default=DEFAULT_BAND, help="band half-width about the observed front [m]")
    ap.add_argument("--tolerance", type=float, default=0.6, help="largest |t_yr - epoch year| scored [yr]")
    ap.add_argument("--out", default=None, help="JSON output")
    args = ap.parse_args()
    from firedrake.petsc import PETSc
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    record = Record(args.kit, comm)
    results = []
    for path in args.states:
        res = score_state(path, record, args.band, args.tolerance, comm)
        if comm.rank == 0:
            results.append(res)
            if res.get("epoch") is None:
                PETSc.Sys.Print(f"  {res['file']}: t={res['t_yr']:.2f}, no Greene epoch within {args.tolerance} yr")
            else:
                f = res["front"]
                PETSc.Sys.Print(f"  {res['file']}: t={res['t_yr']:.2f} vs epoch {res['epoch_year']:.2f}: front rms "
                                f"{f['rms_m'] / 1e3:.2f} km, bias {f['bias_m'] / 1e3:+.2f} km, seaward {f['seaward_share']:.2f}, "
                                f"within 2 km {f['within_2km']:.2f}, front {res['front_km']:.0f} km")
    if comm.rank == 0:
        scored = [r for r in results if r.get("epoch") is not None]
        rms = [r["front"]["rms_m"] for r in scored if r["front"].get("n")]
        summary = {"mean_rms_m": float(np.mean(rms)) if rms else None,
                   "final_rms_m": float(rms[-1]) if rms else None,
                   "mean_bias_m": float(np.mean([r["front"]["bias_m"] for r in scored if r["front"].get("n")])) if rms else None,
                   "epochs": [r["epoch_year"] for r in scored], "band_m": args.band}
        PETSc.Sys.Print(f"summary: mean front rms {summary['mean_rms_m'] and summary['mean_rms_m'] / 1e3:.2f} km over "
                        f"{len(scored)} epochs, final {summary['final_rms_m'] and summary['final_rms_m'] / 1e3:.2f} km")
        if args.out:
            with open(args.out, "w") as f:
                json.dump({"summary": summary, "states": results}, f, indent=1)


if __name__ == "__main__":
    main()
