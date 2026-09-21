#!/usr/bin/env python3
"""Render the open blockers of icepack/ismip7 as the committed NOW.md index.

State and evidence live in the issues. Status and the claim live on the board
(https://github.com/users/dlilien/projects/1). This file is the offline mirror,
for a login node with no `gh`.

    python antarctica/scripts/build_now.py --write
    python antarctica/scripts/build_now.py --check          # exits 1 on drift
    python antarctica/scripts/build_now.py --from-json dump.json --write
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]

REPO = "icepack/ismip7"
BOARD = "https://github.com/users/dlilien/projects/1"
BOARD_OWNER = "dlilien"
BOARD_NUMBER = "1"
BOARD_ID = "PVT_kwHOAHDh9s4BkIgk"

SEVERITY = ("blocks-submission", "owed", "after-deadline")
UNBLOCKER = ("needs-decision", "needs-run", "needs-check", "needs-upstream")
UNVERIFIED = ("unverifiable-nots", "unverifiable-run")
SOURCES = ("src:readiness", "src:submission-readme", "src:matrix-status",
           "src:topic-doc", "src:runbook", "src:open-pr")
# Decisions first: they have the longest lead time and they gate the runs.
UNBLOCKER_ORDER = {name: n for n, name in enumerate(UNBLOCKER)}


def _run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"{' '.join(cmd)} failed: {r.stderr.strip()[:300]}")
    return r.stdout


# `gh project item-list` flattens a multi-select to nothing, so the board is
# read through the API, which returns every option of the Sites field.
BOARD_QUERY = """
query($p:ID!, $after:String){
  node(id:$p){ ... on ProjectV2 { items(first:100, after:$after){
    pageInfo{ hasNextPage endCursor }
    nodes{
      content{ ... on Issue { number } }
      fieldValues(first:20){ nodes{
        __typename
        ... on ProjectV2ItemFieldSingleSelectValue {
          name field{ ... on ProjectV2FieldCommon { name } } }
        ... on ProjectV2ItemFieldMultiSelectValue {
          options{ name } field{ ... on ProjectV2FieldCommon { name } } }
      } }
    } } } }
}
"""


def fetch():
    """Issues from the repository, status and the claim from the board."""
    issues = json.loads(_run([
        "gh", "issue", "list", "--repo", REPO, "--state", "open", "--limit", "200",
        "--json", "number,title,labels,assignees,milestone,updatedAt"]))
    nodes, after = [], None
    while True:
        args = ["gh", "api", "graphql", "-f", f"query={BOARD_QUERY}",
                "-F", f"p={BOARD_ID}"]
        if after:
            args += ["-F", f"after={after}"]
        page = json.loads(_run(args))["data"]["node"]["items"]
        nodes += page["nodes"]
        if not page["pageInfo"]["hasNextPage"]:
            break
        after = page["pageInfo"]["endCursor"]
    return {"issues": issues, "board": nodes}


def index(dump):
    cards = {}
    for it in dump["board"]:
        content = it.get("content") or {}
        if content.get("number") is None:
            continue
        got = {}
        for fv in it["fieldValues"]["nodes"]:
            field = (fv.get("field") or {}).get("name")
            if not field:
                continue
            got[field] = ([o["name"] for o in fv["options"]]
                          if fv.get("options") is not None else fv.get("name"))
        cards[content["number"]] = {
            "status": got.get("Status") or "",
            "owner": got.get("Owner") or "",
            "sites": got.get("Sites") or [],
        }
    rows = []
    for iss in dump["issues"]:
        names = [lab["name"] for lab in iss["labels"]]
        card = cards.get(iss["number"], {})
        rows.append({
            "n": iss["number"],
            "title": iss["title"],
            "labels": names,
            "severity": next((s for s in SEVERITY if s in names), "owed"),
            "unblocker": next((u for u in UNBLOCKER if u in names), ""),
            "sources": [s for s in SOURCES if s in names],
            "unverified": next((u for u in UNVERIFIED if u in names), ""),
            "fresh": "fresh-24h" in names,
            "status": card.get("status", ""),
            "owner": card.get("owner", ""),
            "sites": card.get("sites", []),
            "updated": iss["updatedAt"],
        })
    rows.sort(key=lambda r: (UNBLOCKER_ORDER.get(r["unblocker"], 9), r["n"]))
    return rows


def table(rows):
    out = ["| # | item | unblocked by | owner | sites | state |",
           "|---|---|---|---|---|---|"]
    for r in rows:
        state = r["unverified"] or ("fresh, unverified by design" if r["fresh"]
                                    else "verified open")
        out.append(f'| [{r["n"]}](https://github.com/{REPO}/issues/{r["n"]}) '
                   f'| {r["title"]} | {r["unblocker"] or "n/a"} '
                   f'| {r["owner"] or "unassigned"} '
                   f'| {", ".join(r["sites"]) or "n/a"} | {state} |')
    return out


def render(rows, stamp):
    blocking = [r for r in rows if r["severity"] == "blocks-submission"]
    owed = [r for r in rows if r["severity"] == "owed"]
    later = [r for r in rows if r["severity"] == "after-deadline"]
    claimed = [r for r in rows if r["status"] == "Claimed"]
    unver = [r for r in rows if r["unverified"]]

    L = [
        "# Now",
        "",
        "Generated by `antarctica/scripts/build_now.py` from the open issues of",
        f"`{REPO}`. Do not edit by hand: run `make -C antarctica now`.",
        "",
        f"Status and the claim live on the board, {BOARD}. Evidence lives in the",
        "issues. Reasoning and measured numbers stay in the topic docs. This file",
        "is an index.",
        "",
        "**Claim the card before you start a run or an inversion.** A card in",
        "Claimed for more than 24 hours with no comment reads as unclaimed.",
        "",
        f"Generated: {stamp}. {len(blocking)} blocking, {len(owed)} owed, "
        f"{len(later)} after the deadline, {len(unver)} unverified.",
        "",
        "## Blocking the submission",
        "",
        "Decisions first: they have the longest lead time and they gate the runs.",
        "",
    ]
    L += table(blocking) if blocking else ["None."]
    L += ["", "## Owed", ""]
    L += table(owed) if owed else ["None."]
    L += ["", "## Claimed now", "",
          "Work in flight. Do not duplicate it.", ""]
    if claimed:
        L += ["| # | item | owner | sites | last touched |", "|---|---|---|---|---|"]
        for r in claimed:
            L.append(f'| [{r["n"]}](https://github.com/{REPO}/issues/{r["n"]}) '
                     f'| {r["title"]} | {r["owner"] or "unassigned"} '
                     f'| {", ".join(r["sites"]) or "n/a"} | {r["updated"][:10]} |')
    else:
        L.append("Nothing is claimed.")
    L += ["", "## Unverified", "",
          "Filed open because a check could not settle them from this checkout.",
          "`unverifiable-nots` needs someone at Rice; `unverifiable-run` needs a",
          "production run or a group decision.", ""]
    if unver:
        L += ["| # | item | why |", "|---|---|---|"]
        for r in unver:
            L.append(f'| [{r["n"]}](https://github.com/{REPO}/issues/{r["n"]}) '
                     f'| {r["title"]} | `{r["unverified"]}` |')
    else:
        L.append("None.")
    L += ["", "## Sources", ""]
    counts = {s: sum(1 for r in rows if s in r["sources"]) for s in SOURCES}
    L += ["| source | items |", "|---|---|"]
    for s, c in counts.items():
        L.append(f"| `{s}` | {c} |")
    L += ["", "## After the deadline", ""]
    if later:
        for r in later:
            L.append(f'- [{r["n"]}](https://github.com/{REPO}/issues/{r["n"]}) '
                     f'{r["title"]}')
    else:
        L.append("None.")
    L.append("")
    return "\n".join(L)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=_ROOT / "NOW.md")
    p.add_argument("--from-json", type=Path,
                   help="a saved dump, for a host where gh is not authenticated")
    p.add_argument("--save-json", type=Path,
                   help="write the fetched dump, to carry to such a host")
    p.add_argument("--write", action="store_true")
    p.add_argument("--check", action="store_true",
                   help="exit 1 when the committed file differs from the render")
    a = p.parse_args()

    dump = json.loads(a.from_json.read_text()) if a.from_json else fetch()
    if a.save_json:
        a.save_json.write_text(json.dumps(dump, indent=1))
    rows = index(dump)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    if a.check:
        if not a.output.exists():
            print(f"{a.output} is missing; run make -C antarctica now", file=sys.stderr)
            return 1
        # The stamp is wall clock, so compare everything else.
        def strip(t):
            return "\n".join(l for l in t.splitlines()
                             if not l.startswith("Generated: "))
        if strip(a.output.read_text()) != strip(render(rows, stamp)):
            print(f"{a.output} is stale; run make -C antarctica now", file=sys.stderr)
            return 1
        print(f"{a.output} matches the open issues")
        return 0

    text = render(rows, stamp)
    if a.write:
        a.output.write_text(text)
        print(f"wrote {a.output} from {len(rows)} open issues")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
