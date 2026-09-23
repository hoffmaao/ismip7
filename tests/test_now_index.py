r"""The blocker board is the only place a blocker's status lives.

This repository has no CI, so the unit suite is the one gate a change passes
(`.no-mistakes.yaml`). These tests are what makes "file an issue before you
write the paragraph" hold: a documented open item with no issue reference fails
the suite, and the failure says which command files it.

`NOW.md` is generated from the open issues. The drift check and the reference
check need `gh`, so they skip where it is absent, exactly as the Firedrake
tests skip where Firedrake is absent.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BUILD = REPO / "antarctica" / "scripts" / "build_now.py"
INDEX = REPO / "NOW.md"

# Generated files carry their own provenance and are never hand annotated.
GENERATED = ("NOW.md", "TIMING_MATRIX")
GENERATED_FILES = ("antarctica/reports/core",)

# A line matching one of these states that something is open. It must carry an
# issue reference, on itself or on the heading above it.
MARKERS = (
    (re.compile(r"\bTODO\b"), "TODO"),
    (re.compile(r"\bFIXME\b"), "FIXME"),
    (re.compile(r"\bXXX\b"), "XXX"),
    (re.compile(r"^\s*[-*]\s*\[ \]"), "an unchecked box"),
    (re.compile(r"\[ \]\s+\S"), "an unchecked box"),
    (re.compile(r"\*\*\[confirm(?:, draft)?\]\*\*"), "an unnumbered [confirm]"),
    (re.compile(r"⚠️\s*OPEN"), "an OPEN banner"),
    (re.compile(r"^##+ Open:"), "an Open: heading"),
)
REF = re.compile(r"\(issue #\d+\)|\[confirm(?:, draft)? #\d+\]")
# The legend that explains the checkbox vocabulary is not itself an open item.
# A legend that explains the vocabulary is not itself an open item.
EXEMPT = re.compile(r"`\[x\]`|`\[ \]`|`\[~\]`|item marked|items are marked|marked \*\*\[confirm")


def tracked_markdown():
    out = subprocess.run(["git", "-C", str(REPO), "ls-files", "*.md"],
                         capture_output=True, text=True, check=True).stdout
    for rel in out.split():
        if any(g in rel for g in GENERATED):
            continue
        if any(rel.startswith(f) for f in GENERATED_FILES):
            continue
        yield rel


def offenders():
    found = []
    for rel in tracked_markdown():
        lines = (REPO / rel).read_text().splitlines()
        heading = ""
        for n, line in enumerate(lines, 1):
            if line.startswith("#"):
                heading = line
            if EXEMPT.search(line):
                continue
            for pat, what in MARKERS:
                if pat.search(line) and not REF.search(line) and not REF.search(heading):
                    found.append((rel, n, what, line.strip()[:70]))
                    break
    return found


def test_every_documented_open_item_carries_an_issue():
    bad = offenders()
    if bad:
        report = "\n".join(f"  {r}:{n} carries {w}: {t}" for r, n, w, t in bad)
        pytest.fail(
            "A documented open item has no issue on the blocker board:\n"
            f"{report}\n\n"
            "File it, then put the number on the line:\n"
            "  gh issue create --repo icepack/ismip7 --template blocker.yml\n"
            "  ... then append (issue #NN), or write [confirm #NN]\n"
            "Status and the claim live on the board; this file is an index."
        )


def test_the_detector_would_catch_a_new_marker(tmp_path):
    r"""The guard is worth nothing if it cannot fail."""
    probe = "- [ ] a new blocker nobody filed"
    assert any(p.search(probe) for p, _ in MARKERS)
    assert not REF.search(probe)
    assert REF.search(probe + " (issue #42)")
    # A bare forum-thread number is not a board reference.
    assert not REF.search("departs from the climatology (#48).")


def _gh_ready():
    if shutil.which("gh") is None:
        return False
    return subprocess.run(["gh", "auth", "status"],
                          capture_output=True).returncode == 0


def _board_ready():
    r"""The board is a ProjectV2, which the token reads only with `read:project`."""
    if not _gh_ready():
        return False
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location("build_now", BUILD)
    build_now = module_from_spec(spec)
    spec.loader.exec_module(build_now)
    return subprocess.run(
        ["gh", "api", "graphql", "-F", f"p={build_now.BOARD_ID}", "-f",
         "query=query($p:ID!){node(id:$p){... on ProjectV2{title}}}"],
        capture_output=True).returncode == 0


@pytest.mark.skipif(not _board_ready(),
                    reason="gh cannot read the project board "
                           "(gh auth refresh -s read:project)")
def test_now_index_matches_the_open_issues():
    r = subprocess.run([sys.executable, str(BUILD), "--check"],
                       capture_output=True, text=True)
    # The board is a user project; reading it needs the read:project scope,
    # which a default `gh auth login` token does not carry.
    if r.returncode != 0 and "read:project" in r.stderr:
        pytest.skip("gh token lacks read:project; run gh auth refresh -s read:project")
    assert r.returncode == 0, r.stderr or r.stdout


@pytest.mark.skipif(not _gh_ready(), reason="gh is absent or not authenticated")
def test_every_reference_resolves_to_an_open_issue():
    import json
    open_now = {i["number"] for i in json.loads(subprocess.run(
        ["gh", "issue", "list", "--repo", "icepack/ismip7", "--state", "open",
         "--limit", "200", "--json", "number"],
        capture_output=True, text=True, check=True).stdout)}
    stale = []
    for rel in tracked_markdown():
        for n, line in enumerate((REPO / rel).read_text().splitlines(), 1):
            for num in re.findall(r"\(issue #(\d+)\)|\[confirm(?:, draft)? #(\d+)\]", line):
                num = int(num[0] or num[1])
                if num not in open_now:
                    stale.append(f"  {rel}:{n} points at #{num}, which is closed")
    assert not stale, ("A document still describes something as open whose issue "
                       "is closed:\n" + "\n".join(stale))
