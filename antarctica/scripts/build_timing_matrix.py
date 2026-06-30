#!/usr/bin/env python3
"""
Aggregate per-run timing JSON records into TIMING_MATRIX.md.

Reads results/timing/timing_*.json (one file per mesh/core-count combo)
and writes a markdown table pivoting resolution vs. wall-clock time.

Usage:
    python scripts/build_timing_matrix.py
"""

import glob
import json
import os
from datetime import datetime, timezone

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIMING_DIR = os.path.join(_ROOT, "results", "timing")
OUT_FN = os.path.join(_ROOT, "TIMING_MATRIX.md")


def load_records():
    records = []
    for path in sorted(glob.glob(os.path.join(TIMING_DIR, "timing_*.json"))):
        with open(path) as f:
            records.append(json.load(f))
    return records


def main():
    records = load_records()
    if not records:
        raise SystemExit(
            f"No timing records found in {TIMING_DIR}/. "
            "Run `make timing` (or at least the transient target) first."
        )

    core_counts = sorted({r["ncores"] for r in records})
    rows = sorted(
        {(r["lc"], r["lc_coarse"]) for r in records},
        key=lambda x: (x[0], x[1]),
    )

    # Index: (lc, lc_coarse, ncores) -> record
    index = {(r["lc"], r["lc_coarse"], r["ncores"]): r for r in records}

    # Mesh metadata (same for all core counts at a given resolution)
    meta = {}
    for lc, lcc in rows:
        for nc in core_counts:
            rec = index.get((lc, lcc, nc))
            if rec:
                meta[(lc, lcc)] = (rec["vertices"], rec["cells"])
                break

    header = (
        "| LC (m) | LC_coarse (m) | Vertices | Cells | "
        + " | ".join(f"{nc} cores (s)" for nc in core_counts)
        + " |"
    )
    ncols = 4 + len(core_counts)
    sep = "|" + "|".join(["---"] * ncols) + "|"

    lines = [
        "# Antarctica timing matrix",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "5-year transient runs (`2015`–`2020`, `dt=1.0`), zero SMB/melt forcing, "
        "warm-started from a single LC=2500 inversion (cross-mesh interpolated).",
        "",
        header,
        sep,
    ]

    for lc, lcc in rows:
        verts, cells = meta.get((lc, lcc), ("—", "—"))
        times = []
        for nc in core_counts:
            rec = index.get((lc, lcc, nc))
            times.append(f"{rec['run_seconds']:.1f}" if rec else "—")
        lines.append(
            f"| {lc} | {lcc} | {verts} | {cells} | "
            + " | ".join(times)
            + " |"
        )

    per_step_header = (
        "| LC (m) | LC_coarse (m) | Vertices | Cells | "
        + " | ".join(f"{nc} cores (s/step)" for nc in core_counts)
        + " |"
    )

    lines.extend(
        [
            "",
            "## Per-step timing",
            "",
            per_step_header,
            sep,
        ]
    )

    for lc, lcc in rows:
        verts, cells = meta.get((lc, lcc), ("—", "—"))
        times = []
        for nc in core_counts:
            rec = index.get((lc, lcc, nc))
            times.append(f"{rec['seconds_per_step']:.2f}" if rec else "—")
        lines.append(
            f"| {lc} | {lcc} | {verts} | {cells} | "
            + " | ".join(times)
            + " |"
        )

    lines.append("")
    with open(OUT_FN, "w") as f:
        f.write("\n".join(lines))
    print(f"Wrote {OUT_FN} ({len(rows)} resolutions x {len(core_counts)} core counts)")


if __name__ == "__main__":
    main()
