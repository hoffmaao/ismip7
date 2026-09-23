#!/usr/bin/env python3
r"""Render the tracked run log from one JSON record per simulation.

Every simulation this project runs -- an inversion, a historical, a control, a
projection, the OCX experiment, a calibration, a forward test, an output
rehearsal -- leaves a record in ``antarctica/runlog/<id>.json``, and this
renders them into ``antarctica/reports/SIMULATIONS.md``. One file per run so
two sites never conflict in the same file, and JSON so the fields can be
checked rather than eyeballed.

Why here rather than in a spreadsheet: the results themselves are gitignored
(checkpoints are hundreds of MB, timeseries are regenerated), so the record is
the only reviewable trace of a run that reaches the repository, exactly as the
per-core reports of ``core_report.py`` are. A record is reviewed in a pull
request beside the code the run used, and ``--check`` refuses a stale render.

``core_report.py`` stays what it is: the deep per-core report of one core
experiment, with its budget rows and audit. This is the index across all of
them, and across the runs that are not core experiments.

The group also keeps a shared progress sheet, which is a second reader of the
same records rather than a second place to type them: ``--csv`` writes them
flat for it to import. A hand-maintained copy of this information is exactly
what this file exists to replace, so the sheet reads the repository.

Usage:
    python antarctica/scripts/build_runlog.py --write
    python antarctica/scripts/build_runlog.py --check       # exits 1 on drift
    python antarctica/scripts/build_runlog.py --csv antarctica/results/runlog.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent
_ANT = _SCRIPTS.parent
RUNLOG_DIR = _ANT / "runlog"
OUTPUT = _ANT / "reports" / "SIMULATIONS.md"

#: Required of every record. Everything else is optional and rendered when set,
#: so a planned run is one short file and a finished one carries its whole
#: provenance.
REQUIRED = ("id", "task", "title", "status", "institution")

#: Task kinds, in the order the report groups them: what a run has to be
#: before it can be the next thing.
TASKS = ("inversion", "calibration", "test", "historical", "control",
         "projection", "ocx", "output")

#: Recognised statuses. ``superseded`` is kept because a result that has been
#: invalidated must stay in the log saying so, not disappear from it.
STATUSES = ("planned", "queued", "running", "stopped", "done", "superseded")

#: Every field, in the order a record's detail block lists them, with the
#: heading each one carries.
FIELDS = (
    ("id", "Record"),
    ("institution", "Institution"),
    ("title", "Simulation"),
    ("status", "Status"),
    ("task", "Task type"),
    ("exp_id", "ISMIP7 exp id"),
    ("esm", "ESM"),
    ("scenario", "Scenario"),
    ("period", "Period (yr)"),
    ("friction", "Friction law"),
    ("mesh", "Mesh"),
    ("initial_state", "Initial state / MAP"),
    ("branch_from", "Branch from"),
    ("forcing", "Forcing versions"),
    ("melt", "Melt: K, slope, deltaT"),
    ("front", "Calving front, collapse"),
    ("apparent_mb", "Apparent MB"),
    ("dt", "dt (yr)"),
    ("site", "Site / partition"),
    ("resources", "Ranks / memory"),
    ("jobs", "Job ids"),
    ("code", "Code"),
    ("started", "Started"),
    ("finished", "Finished"),
    ("cost", "Cost per model year"),
    ("results", "Results path"),
    ("audit", "Audit"),
    ("ismip7_output", "ISMIP7 output written"),
    ("checker", "Regridded, isschecker"),
    ("scalars", "Scalars processed"),
    ("uploaded", "Uploaded"),
    ("notes", "Notes"),
)

#: The summary table's columns. A record's own detail block carries the rest;
#: a table of every field is unreadable at any width.
SUMMARY = (("title", "Simulation"), ("status", "Status"), ("mesh", "Mesh"),
           ("site", "Site"), ("started", "Started"), ("finished", "Finished"),
           ("audit", "Headline result"))


def load(directory=RUNLOG_DIR):
    r"""Every record in ``directory``, validated, sorted by task then id."""
    records, seen = [], {}
    for path in sorted(directory.glob("*.json")):
        try:
            record = json.loads(path.read_text())
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path.name}: not valid JSON: {exc}") from exc
        if not isinstance(record, dict):
            raise ValueError(f"{path.name}: not a JSON object")
        missing = [k for k in REQUIRED if not record.get(k)]
        if missing:
            raise ValueError(f"{path.name}: missing {', '.join(missing)}")
        if record["task"] not in TASKS:
            raise ValueError(
                f"{path.name}: task {record['task']!r} is not one of "
                f"{', '.join(TASKS)}")
        if record["status"] not in STATUSES:
            raise ValueError(
                f"{path.name}: status {record['status']!r} is not one of "
                f"{', '.join(STATUSES)}")
        known = {k for k, _ in FIELDS}
        unknown = sorted(set(record) - known)
        if unknown:
            raise ValueError(
                f"{path.name}: unknown field(s) {', '.join(unknown)}; add them "
                f"to FIELDS in build_runlog.py if they belong in the log")
        if record["id"] != path.stem:
            raise ValueError(
                f"{path.name}: id {record['id']!r} does not match the filename")
        if record["id"] in seen:
            raise ValueError(f"duplicate id {record['id']!r}")
        seen[record["id"]] = path
        records.append(record)
    records.sort(key=lambda r: (TASKS.index(r["task"]), r["id"]))
    return records


def _cell(value):
    r"""One field as table text: a list joins, an empty value is a dash, and a
    pipe would end the column early."""
    if value is None or value == "" or value == []:
        return "-"
    if isinstance(value, (list, tuple)):
        value = " ".join(str(v) for v in value)
    return str(value).replace("|", "\\|").replace("\n", " ")


def render(records):
    out = [
        "# ISMIP7 Antarctica simulations",
        "",
        "Generated by `antarctica/scripts/build_runlog.py` from the records in",
        "`antarctica/runlog/`. Do not edit by hand: add or update a record and",
        "run `make -C antarctica runlog`.",
        "",
        "One record per simulation, of any kind. The results themselves are",
        "gitignored, so these records and the per-core reports beside them are",
        "the trace a run leaves in the repository. A core experiment also gets",
        "its full report from `core_report.py`; this is the index.",
        "",
        f"{len(records)} records.",
        "",
    ]
    counts = {}
    for r in records:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    if counts:
        out += ["Status: " + ", ".join(
            f"{counts[s]} {s}" for s in STATUSES if s in counts) + ".", ""]

    for task in TASKS:
        group = [r for r in records if r["task"] == task]
        if not group:
            continue
        out += [f"## {task.capitalize()}", ""]
        out += ["| " + " | ".join(h for _, h in SUMMARY) + " |",
                "|" + "|".join("---" for _ in SUMMARY) + "|"]
        for r in group:
            out.append("| " + " | ".join(
                _cell(r.get(k)) for k, _ in SUMMARY) + " |")
        out.append("")

    out += ["## Detail", ""]
    for r in records:
        out += [f"### {r['id']}", "", f"{r['title']} ({r['status']}), "
                f"{r['institution']}.", ""]
        for key, heading in FIELDS:
            if key in ("id", "title", "status", "institution"):
                continue
            value = r.get(key)
            if value in (None, "", []):
                continue
            out.append(f"- **{heading}:** {_cell(value)}")
        out.append("")
    return "\n".join(out).rstrip() + "\n"


def _csv_cell(value):
    r"""One field as sheet text: a list joins, an empty value stays empty, and
    everything else is the record's own text, which the csv module quotes."""
    if value is None or value == [] or value == "":
        return ""
    if isinstance(value, (list, tuple)):
        return " ".join(str(v) for v in value)
    return str(value)


def write_csv(records, path):
    r"""The records flat, one row each, for the group's shared progress sheet.

    Every field of :data:`FIELDS` in that order, so the sheet's header is the
    generator's own list and a field added there reaches the sheet on the next
    import. Deliberately outside ``--check``: the markdown is the tracked
    render, and this is a copy made on demand for a reader outside the
    repository, so that reader is never a second place to type any of it."""
    import csv

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([heading for _, heading in FIELDS])
        for record in records:
            writer.writerow([
                _csv_cell(record.get(key)) for key, _ in FIELDS])
    return len(records)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--runlog", type=Path, default=RUNLOG_DIR)
    p.add_argument("--output", type=Path, default=OUTPUT)
    p.add_argument("--write", action="store_true")
    p.add_argument("--check", action="store_true",
                   help="exit 1 when the committed file differs from the render")
    p.add_argument("--csv", type=Path,
                   help="also write the records flat, for the shared sheet")
    a = p.parse_args(argv)

    try:
        records = load(a.runlog)
    except ValueError as exc:
        print(f"runlog: {exc}", file=sys.stderr)
        return 2
    text = render(records)

    if a.csv:
        write_csv(records, a.csv)
        print(f"wrote {a.csv} with {len(records)} records")

    if a.check:
        if not a.output.exists():
            print(f"{a.output} is missing; run make -C antarctica runlog",
                  file=sys.stderr)
            return 1
        if a.output.read_text() != text:
            print(f"{a.output} is stale; run make -C antarctica runlog",
                  file=sys.stderr)
            return 1
        print(f"{a.output} matches the {len(records)} records")
        return 0

    if a.write:
        a.output.write_text(text)
        print(f"wrote {a.output} from {len(records)} records")
    elif not a.csv:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
