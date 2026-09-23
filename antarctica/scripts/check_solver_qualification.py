#!/usr/bin/env python3
r"""Fail-loud gate for the transient solver qualification probes.

The timing driver records requested and completed steps separately because an
exhausted rescue ladder saves and returns cleanly.  This checker requires the
whole interval, the intended number of CSV rows, no diverged diagnostic solve,
and a mass-budget residual that rounds to zero in the persisted CSV.
"""

import argparse
import csv
import json
import math


def _is_diverged(reason):
    if reason.startswith("DIVERGED"):
        return True
    try:
        return int(reason) < 0
    except ValueError:
        return False


def qualification_verdict(record, rows, min_steps, residual_tol_gt=5e-5,
                          kinds=("qualification",)):
    """``(ok, detail)`` for one lane record and its timeseries rows.

    The gate `make qualify` runs, as a function so the map-check manager can
    hold its own 10-step lanes to the same rule: the whole interval, the
    intended number of rows, no diverged solve, a residual that rounds to
    zero. ``kinds`` names the record kinds accepted.
    """
    if record.get("timing_kind") not in kinds:
        return False, (
            f"not a {'/'.join(kinds)} record "
            f"(timing_kind={record.get('timing_kind')!r})"
        )
    if record.get("nsteps", 0) < min_steps:
        return False, (
            f"record requested {record.get('nsteps', 0)} steps; "
            f"need at least {min_steps}"
        )
    if record.get("completed_steps", 0) < min_steps:
        return False, (
            f"record completed {record.get('completed_steps', 0)} steps; "
            f"need at least {min_steps}"
        )

    summary = record.get("diagnostic_solve_summary", {})
    reasons = summary.get("reason_counts", {})
    diverged = {name: count for name, count in reasons.items()
                if _is_diverged(name)}
    if diverged:
        return False, f"diagnostic solve divergence recorded: {diverged}"

    transport_summary = record.get("transport_solve_summary", {})
    transport_reasons = transport_summary.get("reason_counts", {})
    transport_diverged = {
        name: count for name, count in transport_reasons.items()
        if _is_diverged(name)
    }
    if transport_diverged:
        return False, f"transport KSP divergence recorded: {transport_diverged}"
    if transport_summary.get("count", 0) < min_steps:
        return False, (
            "record contains only "
            f"{transport_summary.get('count', 0)} transport solves; "
            f"need at least {min_steps}"
        )
    transport_mass_resid = transport_summary.get("mass_residual_gt_max")
    if (
        transport_mass_resid is None
        or not math.isfinite(transport_mass_resid)
        or transport_mass_resid > residual_tol_gt
    ):
        return False, (
            "max |transport mass residual| is "
            f"{transport_mass_resid!r} Gt; limit is {residual_tol_gt:.6g} Gt"
        )

    if len(rows) < min_steps:
        return False, (
            f"timeseries contains {len(rows)} completed steps; "
            f"need at least {min_steps}"
        )

    last_year = float(rows[-1]["year"])
    expected_year = float(record["t_end"])
    if not math.isclose(last_year, expected_year, abs_tol=0.51 * record["dt"]):
        return False, f"timeseries stops at {last_year}, expected {expected_year}"

    residuals = [abs(float(row["resid_gt"])) for row in rows]
    max_residual = max(residuals, default=math.inf)
    if not math.isfinite(max_residual) or max_residual > residual_tol_gt:
        return False, (
            f"max |mass residual| is {max_residual:.6g} Gt; "
            f"limit is {residual_tol_gt:.6g} Gt"
        )

    mode = record.get("diagnostic_solver_mode", "unknown")
    return True, (
        f"{mode}: {len(rows)} transient steps, "
        f"{summary.get('count', 0)} diagnostic solves, "
        f"max |transport/step mass residual|="
        f"{transport_mass_resid:.3g}/{max_residual:.3g} Gt"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--min-steps", type=int, required=True)
    # CSV stores four decimal places, so 5e-5 accepts only rows persisted as
    # 0.0000 (including signed zero), not a small but nonzero rounded residual.
    parser.add_argument("--residual-tol-gt", type=float, default=5e-5)
    parser.add_argument(
        "--kind", action="append", default=None,
        help="record kind(s) accepted; repeatable (default: qualification)",
    )
    args = parser.parse_args()

    with open(args.record) as stream:
        record = json.load(stream)
    with open(args.csv, newline="") as stream:
        rows = list(csv.DictReader(stream))
    ok, detail = qualification_verdict(
        record, rows, args.min_steps, args.residual_tol_gt,
        kinds=tuple(args.kind or ("qualification",)),
    )
    if not ok:
        raise SystemExit(f"FAIL: {detail}")
    print(f"PASS: {detail}")


if __name__ == "__main__":
    main()
