#!/usr/bin/env python
r"""What the gridded melt and calving fields carry and leave out over one
run's annual files, and whether the booking keeps issue 136's rules. Read
only, serial.

    python antarctica/scripts/check_melt_booking.py STEM [--start FILE]
        [--overlap CACHE] [--timeseries CSV] [--first YEAR] [--last YEAR]
        [--tol GT_PER_YR] [--csv OUT] [--report-only]

STEM names the series as the writer does, ``<results>/<run>_ismip7_annual.h5``
(the years are ``<stem>_<year>.h5``). A year's start state is the year before
it in the series; for the first year it is ``--start``, the annual file whose
year ends where the series begins (a projection's historical's last year),
and without it the first year is read as a start state only. ``--overlap`` is
the writer's cache, ``<STEM>.overlap.npz`` by default, and is built in memory
when absent; the writer builds it in the cell order of the series' first year,
and every file is matched to that order by centroid, so a start file written
on another rank count reads right. ``--timeseries`` is the forward's own
``<run>_timeseries.csv``, whose ``calv_gt`` and ``outflux_gtyr`` columns
``licalvf`` has to carry.

Per year, in Gt/yr over map-plane area, ice at 917 kg m-3:

* ``lifmassbf``, the melt of cells with no floating ice at year end, split by
  the cells that booked it: cells holding no ice at either end (inflow that
  melted on arrival), shelf ice afloat when the year began (gone or grounded
  by its end), and the rest (cells that began the year grounded or empty);
* ``left out``, the melt in 8 km pixels with no floating ice at year end
  after the writer's near-flotation rule, which the ``no_floating_ice`` fill
  drops from ``libmassbffl``, and ``refreezing left out``, the refreezing
  there, which ``lifmassbf``'s range excludes (reported only);
* ``reference in sinks``, on cells holding no ice at either end, the positive
  apparent-MB reference that the booked melt and negative SMB could have
  taken, ``min(corr+, M + acabf-)`` (``ismip7_output.net_reference``);
* ``snowfall in sinks``, on the same cells, the positive SMB that the booked
  melt could have taken, ``min(acabf+, M)`` (``ismip7_output.net_snowfall``);
* ``negative reference``, the negative reference booked on cells holding no
  ice at either end, a sink on real ice: the part of their inflow it removed,
  which stays in ``acabf_correction`` and reaches no submitted field;
  ``reference on melting cells``, the reference booked on every cell that
  books melt; and ``reference applied``, the reference booked on every cell,
  which no submitted field carries, so the submitted fields' mass budget
  misses it by that much (all three reported only, for issue #104);
* ``land melt``, the melt booked on land cells (a bed at or above sea level)
  holding no ice at either end;
* ``land ligroundf``, the grounding-line flux booked into those land cells;
* ``licalvf``, and ``licalvf against the timeseries``, its difference from
  the year's ``calv`` plus ``outflux`` in the forward's timeseries (with
  ``--timeseries``);
* ``budget residual``, the year's change in thickness less every flux booked,
  ``acabf``, ``acabf_correction``, ``libmassbffl``, ``lifmassbf`` and
  ``licalvf``: the transport moves ice between cells and nothing else, so it
  is zero when every loss has a field;
* ``split``, the largest difference between the ``lifmassbf`` written and
  ``ismip7_output.split_melt`` applied to the file's own melt and written
  floating mask, in m/yr.

Exit status 0 when every year keeps the left-out melt, the reference and the
snowfall in sinks, the land melt and the land ligroundf at or below ``--tol``,
the budget residual and ``licalvf`` against the timeseries at or below 0.05
Gt/yr, and the split within 1e-9 m/yr; 1 when one does not
(``--report-only`` reports and exits 0); 2 when an input is missing. A series
booked before issue 136's changes of 28 September 2026 fails it: that is what
the gate is for.
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

GT = 917.0 / 1e12                     # m3 of ice to Gt
SPLIT_TOL = 1e-9                      # m/yr
# Gt/yr. The forward holds each advance's transport to a mass residual of
# ISMIP7_MASS_RESIDUAL_TOL_GT, 5e-5 Gt by default, and a year holds 40 steps
# of up to 16 subcycled advances; the timeseries carries four decimals.
BUDGET_TOL = 0.05
MARKERS = (2100, 2200, 2300)
FIELDS = ("libmassbffl", "lifmassbf", "acabf", "acabf_correction", "dlithkdt", "ligroundf",
          "licalvf", "sftgif", "sftflf", "sftgrf", "orog", "lithk", "topg")
GATED = ("left out", "reference in sinks", "snowfall in sinks", "land melt", "land ligroundf")
BUDGET_GATED = ("budget residual", "licalvf against the timeseries")


def year_booking(prev, cur, W):
    r"""One year's quantities (the module docstring's list, less the
    timeseries) from the start and end states' fields, dicts of per-cell
    arrays in the cell order of ``W``, the overlap operator (pixels by cells,
    m2). Pure numpy."""
    import write_ismip7_output as wio
    from icepack2_tools.ismip7_output import split_melt
    area = np.asarray(W.sum(axis=0)).ravel()
    L, F, S, R, C = (cur[k] for k in ("libmassbffl", "lifmassbf", "acabf",
                                      "acabf_correction", "licalvf"))
    B = L + F                                            # the melt before the split
    M = np.maximum(-B, 0.0)
    ice0, ice1 = prev["sftgif"] > 0.5, cur["sftgif"] > 0.5
    none = ~ice0 & ~ice1
    land = cur["topg"] >= 0.0
    afloat0 = ice0 & (prev["sftflf"] > 0.5)
    masks = {k: cur[k].copy() for k in ("sftgif", "sftflf", "sftgrf", "orog", "lithk", "topg")}
    wio.ground_near_flotation(masks)
    floating = masks["sftflf"] > 0.5
    _, lif = split_melt(B, floating)
    left = (W @ floating.astype(float)) <= 0.0
    in_sinks = np.where(none, np.minimum(np.maximum(R, 0.0), M + np.maximum(-S, 0.0)), 0.0)
    snow = np.where(none, np.minimum(np.maximum(S, 0.0), M), 0.0)

    def gt(v, where=True):
        return float((v * area * where).sum()) * GT

    def left_out(v):
        return float((W @ v)[left].sum()) * GT

    return {
        "lifmassbf": -gt(F),
        "lifmassbf: no ice at either end": -gt(F, none),
        "lifmassbf: afloat at start": -gt(F, afloat0),
        "lifmassbf: the rest": -gt(F, ~afloat0 & ~none),
        "left out": left_out(np.maximum(-L, 0.0)),
        "refreezing left out": left_out(np.maximum(L, 0.0)),
        "reference in sinks": gt(in_sinks),
        "snowfall in sinks": gt(snow),
        "negative reference": gt(np.minimum(R, 0.0), none),
        "reference on melting cells": gt(R, M > 0.0),
        "reference applied": gt(R),
        "land melt": gt(M, none & land),
        "land ligroundf": gt(cur["ligroundf"], none & land),
        "licalvf": -gt(C),
        "budget residual": gt(cur["dlithkdt"] - S - R - B - C),
        "split": float(np.max(np.abs(lif - F), initial=0.0)),
    }


def timeseries_losses(path):
    r"""``{year: calv + outflux}`` in Gt over each model year, from the
    forward's timeseries CSV: a row closes the step that ends at its
    ``year``, and a step ending at Y + 1 belongs to year Y. A row repeated by
    a resumed link counts once, as the resume wrote it."""
    with open(path) as f:
        rows = list(csv.DictReader(f))
    by_t = {}
    for r in rows:
        by_t[round(float(r["year"]), 9)] = r                  # the resume's row wins
    t = np.array(sorted(by_t))
    dt = np.diff(t, prepend=t[0] - (t[1] - t[0]) if len(t) > 1 else t[0] - 1.0)
    loss = np.array([float(by_t[k]["calv_gt"]) for k in t])
    loss += np.array([float(by_t[k]["outflux_gtyr"]) for k in t]) * dt
    year = np.ceil(t - 1e-9).astype(int) - 1
    return {int(y): float(loss[year == y].sum()) for y in np.unique(year)}


def failures(rows, tol):
    r"""The gate: every ``(year, quantity, value)`` past its bound."""
    out = [(r["year"], k, r[k]) for r in rows for k in GATED if abs(r[k]) > tol]
    # a year the timeseries lacks compares as NaN, and fails
    out += [(r["year"], k, r[k]) for r in rows for k in BUDGET_GATED
            if k in r and not abs(r[k]) <= BUDGET_TOL]
    out += [(r["year"], "split", r["split"]) for r in rows if r["split"] > SPLIT_TOL]
    return out


def summary(rows):
    r"""Each quantity at its largest and at the marker years the series holds."""
    lines = []
    years = [r["year"] for r in rows]
    shown = [y for y in MARKERS if y in years]
    if years and years[-1] not in shown:
        shown.append(years[-1])
    by_year = {r["year"]: r for r in rows}
    for k in rows[0]:
        if k == "year":
            continue
        top = max(rows, key=lambda r: abs(r[k]))
        at = ", ".join(f"{by_year[y][k]:.4g} in {y}" for y in shown)
        lines.append(f"  {k}: largest {top[k]:.4g} ({top['year']}); {at}")
    return lines


# ---- reading (Firedrake) ---------------------------------------------------

def _load(path, names):
    import firedrake as fd
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        f = {n: chk.load_function(mesh, n).dat.data_ro.copy() for n in names}
        xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(
            fd.SpatialCoordinate(mesh)).dat.data_ro.copy()
    return f, xy, mesh


class _Aligner:
    r"""Every file's cells in one reference order, matched by centroid."""

    def __init__(self, xy):
        self.xy = xy
        self.inv = np.argsort(np.lexsort((xy[:, 1], xy[:, 0])))

    def __call__(self, f, xy):
        perm = np.lexsort((xy[:, 1], xy[:, 0]))[self.inv]
        if xy.shape != self.xy.shape or not np.allclose(xy[perm], self.xy):
            raise ValueError("the files do not share one mesh")
        return {k: v[perm] for k, v in f.items()}


def _years(stem):
    pat = re.compile(re.escape(stem) + r"_(\d+)\.h5$")
    return sorted(int(m.group(1)) for p in glob.glob(glob.escape(stem) + "_*.h5")
                  for m in [pat.match(p)] if m)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stem")
    ap.add_argument("--start")
    ap.add_argument("--overlap")
    ap.add_argument("--timeseries")
    ap.add_argument("--first", type=int)
    ap.add_argument("--last", type=int)
    ap.add_argument("--tol", type=float, default=1e-6, help="Gt/yr (default 1e-6)")
    ap.add_argument("--csv")
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args(argv)

    stem = a.stem[:-3] if a.stem.endswith(".h5") else a.stem
    years = _years(stem)
    cache = a.overlap or f"{stem}.h5.overlap.npz"
    missing = [p for p in (a.start, a.timeseries) if p and not os.path.exists(p)]
    if not years or missing:
        print(f"no annual files under {stem}_<year>.h5" if not years
              else f"missing: {', '.join(missing)}", file=sys.stderr)
        return 2
    import write_ismip7_output as wio
    first, xy0, mesh0 = _load(f"{stem}_{years[0]}.h5", FIELDS)
    align = _Aligner(xy0)
    W = wio.overlap_operator(mesh0 if not os.path.exists(cache) else None,
                             cache if os.path.exists(cache) else None)
    if W.shape[1] != len(xy0):
        print(f"{cache} covers {W.shape[1]} cells, the series {len(xy0)}", file=sys.stderr)
        return 2
    losses = timeseries_losses(a.timeseries) if a.timeseries else None
    lo = a.first if a.first is not None else years[0]
    hi = a.last if a.last is not None else years[-1]
    todo = [y for y in years if lo <= y <= hi]
    if todo and todo[0] - 1 in years:
        prev = align(*_load(f"{stem}_{todo[0] - 1}.h5", FIELDS)[:2])
    elif a.start:
        prev = align(*_load(a.start, FIELDS)[:2])
    else:
        prev, todo = align(first, xy0), todo[1:]
    rows = []
    for yr in todo:
        cur = align(*_load(f"{stem}_{yr}.h5", FIELDS)[:2])
        r = {"year": yr, **year_booking(prev, cur, W)}
        if losses is not None:
            r["licalvf against the timeseries"] = r["licalvf"] - losses.get(yr, np.nan)
        rows.append(r)
        ts = (f", {r['licalvf against the timeseries']:+.3g} against the timeseries"
              if losses is not None else "")
        print(f"{yr}: lifmassbf {r['lifmassbf']:,.1f} (no ice at either end "
              f"{r['lifmassbf: no ice at either end']:,.1f}, afloat at start "
              f"{r['lifmassbf: afloat at start']:,.1f}, the rest {r['lifmassbf: the rest']:,.1f}), "
              f"left out {r['left out']:.3g}, reference in sinks {r['reference in sinks']:.3g}, "
              f"snowfall in sinks {r['snowfall in sinks']:.3g}, negative reference "
              f"{r['negative reference']:.3g}, land melt {r['land melt']:.3g}, land ligroundf "
              f"{r['land ligroundf']:.3g}, licalvf {r['licalvf']:,.1f}{ts}, budget residual "
              f"{r['budget residual']:+.3g} Gt/yr; split {r['split']:.1e} m/yr", flush=True)
        prev = cur
    if rows:
        print(f"{os.path.basename(stem)}, {rows[0]['year']} to {rows[-1]['year']}:")
        for line in summary(rows):
            print(line)
    if a.csv and rows:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    bad = failures(rows, a.tol)
    for where, what, value in bad:
        print(f"FAIL {where}: {what} {value:.4g}")
    print(f"check_melt_booking: {len(bad)} failures "
          f"({'reported only' if a.report_only else 'gated'})")
    return 1 if bad and not a.report_only else 0


if __name__ == "__main__":
    sys.exit(main())
