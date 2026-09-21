#!/usr/bin/env python3
r"""Is the bundled ISMIP7 variable request still the checker's?

    python antarctica/scripts/audit_variable_request.py [--ref main]

``icepack2_tools/ismip7_variable_request.csv`` is a copy of the table the
compliance checker reads (``isschecker/data/ISMIP7_variable_request.csv`` in
``ismip/ISM_SimulationChecker``), and the writer takes its units, standard
names, FL/ST types and fill policies from it. Upstream edits it as the forum
finds problems: the bounds were widened in June and July (discussions #16,
#22, #23), the file moved out of ``conventions/`` in August (#5), and
``range_severity`` changed on 17 September (#46). A copy that falls behind
writes files the current checker reads differently, so this prints what has
moved since the copy was taken, row by row, and exits 1 if anything has.

``ismip7_variable_request.source.json`` beside the copy records the tag and
commit it came from. To take a new copy, download the file at the new tag
over the old one and update that record in the same commit.

Columns the writer reads are marked ``WRITER``: a change there changes the
submission files. The rest (the bounds and severities) only change what the
checker says about them.
"""
import argparse
import csv
import io
import json
import os
import sys
import urllib.request

_TOOLS = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "icepack2_tools")
VENDORED = os.path.join(_TOOLS, "ismip7_variable_request.csv")
SOURCE = os.path.join(_TOOLS, "ismip7_variable_request.source.json")
RAW = "https://raw.githubusercontent.com/{repository}/{ref}/{path}"
KEY = "Variable Name"
# what write_ismip7_output.py and ismip7_output.py take from the table
WRITER_COLUMNS = ("Type", "units", "long_name", "standard_name", "fill_policy", "Dim")


def rows(text):
    return {r[KEY]: r for r in csv.DictReader(io.StringIO(text))}


def differences(ours, theirs):
    r"""``[(variable, column, ours, theirs)]``; a variable on one side only
    has column ``None``."""
    out = []
    for var in sorted(set(ours) | set(theirs)):
        if var not in ours or var not in theirs:
            out.append((var, None, "present" if var in ours else "absent",
                        "present" if var in theirs else "absent"))
            continue
        for col in sorted(set(ours[var]) | set(theirs[var])):
            a, b = ours[var].get(col), theirs[var].get(col)
            if a != b:
                out.append((var, col, a, b))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default="main", help="upstream branch, tag or commit (default %(default)s)")
    a = ap.parse_args()
    with open(SOURCE) as f:
        src = json.load(f)
    with open(VENDORED, newline="") as f:
        ours = rows(f.read())
    url = RAW.format(repository=src["repository"], ref=a.ref, path=src["path"])
    with urllib.request.urlopen(url, timeout=60) as r:
        theirs = rows(r.read().decode())
    print(f"bundled: {src['repository']} {src['tag']} ({src['commit'][:10]}, {src['committed']}), "
          f"{len(ours)} variables\nupstream: {a.ref}, {len(theirs)} variables")
    diffs = differences(ours, theirs)
    for var, col, mine, upstream in diffs:
        if col is None:
            print(f"  {var}: {mine} here, {upstream} upstream")
        else:
            mark = "  WRITER" if col in WRITER_COLUMNS else ""
            print(f"  {var}.{col}: {mine!r} here, {upstream!r} upstream{mark}")
    print("current" if not diffs else f"{len(diffs)} difference(s): take a new copy")
    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main())
