r"""The core matrix's wall-retry policy, per experiment kind.

``run_core_matrix.sh`` relaunches a run that stopped short of its target from
its own saved state, up to MAX_ATTEMPTS, passing ``ISMIP7_RESTART`` on every
attempt after the first. Core 11 (OCX) used to be exempt because its driver
cold-started unconditionally; it now reads ``ISMIP7_RESTART`` and auto-resumes
like every other core, so a wall stop must be retried rather than reported as
unretryable.

``archive_stale`` moves superseded output aside before a re-run. It has to
take the ISMIP7 annual series with the rest: a cold start into a populated
series is a hard refusal (``AnnualOutput`` will not overwrite a banked
submission), so a series left behind makes the re-run impossible.

A core that KEEPS its output (the reusable path, where archive_stale never
runs) must resume from its own state on the first attempt too, for the same
reason: cold-starting a partially complete projection would rewrite the years
it already banked.

The script resolves its repository root from ``BASH_SOURCE``, so it runs here
against a sandbox laid out like the repo: the shipped script is copied in
unmodified, ``mpiexec`` is a PATH stub, and the driver is a stub that writes
the timeseries row and the checkpoint attribute the matrix reads back. Every
decision under test is the shipped script's.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MATRIX = REPO / "antarctica" / "scripts" / "run_core_matrix.sh"

pytest.importorskip("h5py")

# Advances to FAKE_STOP_YEAR on a cold start and to FAKE_TARGET once the matrix
# hands it ISMIP7_RESTART, which is what a wall-stopped core does on retry.
DRIVER = r'''
import os, sys
import h5py

# the matrix cd's to its repo root before launching the driver
R = os.path.join(os.getcwd(), "antarctica", "results")
os.makedirs(R, exist_ok=True)
stem = os.environ["FAKE_STEM"]
csv = os.path.join(R, stem + "_timeseries.csv")
h5 = os.path.join(R, stem + "_final.h5")

restart = os.environ.get("ISMIP7_RESTART")
with open(os.path.join(R, "attempts.log"), "a") as f:
    f.write(("restart" if restart else "cold") + "\n")

end = float(os.environ["FAKE_TARGET"] if restart
            else os.environ["FAKE_STOP_YEAR"])
start = float(os.environ["FAKE_START_YEAR"])
new = not os.path.exists(csv)
with open(csv, "a") as f:
    if new:
        f.write("year,vaf\n")
    y = start
    while y <= end + 1e-9:
        f.write(f"{y:.1f},0.0\n")
        y += 1.0
with h5py.File(h5, "w") as h:
    h["/"].attrs["t_yr"] = end
    h["/"].attrs["stalled"] = 0
print(f"Saved: {h5}")
'''

MPIEXEC = '#!/bin/bash\nexec "$FAKE_PYTHON" "$FAKE_DRIVER"\n'


@pytest.fixture
def sandbox(tmp_path):
    r"""A repo-shaped tree holding the shipped matrix script and the stubs."""
    scripts = tmp_path / "antarctica" / "scripts"
    (scripts / "projections").mkdir(parents=True)
    (tmp_path / "antarctica" / "results").mkdir(parents=True)
    (tmp_path / "icepack2_tools").mkdir()
    # PROV_REF: results must postdate it, and the stub writes them after this
    (tmp_path / "icepack2_tools" / "forcing.py").write_text("")

    shutil.copy(MATRIX, scripts / "run_core_matrix.sh")
    (scripts / "projections" / "ocx.py").write_text(DRIVER)

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "mpiexec").write_text(MPIEXEC)
    (bin_dir / "mpiexec").chmod(0o755)
    venv_bin = tmp_path / "venv" / "bin"
    venv_bin.mkdir(parents=True)
    (venv_bin / "python").symlink_to(sys.executable)
    return tmp_path


def _run(sandbox, **env):
    e = dict(os.environ)
    e.update({
        "PATH": f"{sandbox / 'bin'}:{e['PATH']}",
        "VENV": str(sandbox / "venv"),
        "FAKE_PYTHON": sys.executable,
        "FAKE_DRIVER": str(sandbox / "antarctica" / "scripts" / "projections" / "ocx.py"),
        "MAX_LOAD": "100000",          # never gate on load in a test
        "NRANKS": "1",
        "CORES": "11",
        "FRESH": "1",
        "ISMIP7_LC": "32000",
    })
    e.pop("ISMIP7_RUN_TAG", None)
    e.update(env)
    p = subprocess.run(
        ["bash", str(sandbox / "antarctica" / "scripts" / "run_core_matrix.sh")],
        env=e, cwd=str(sandbox), capture_output=True, text=True, timeout=300,
    )
    return p.stdout + p.stderr


def _attempts(sandbox):
    log = sandbox / "antarctica" / "results" / "attempts.log"
    return log.read_text().split() if log.exists() else []


def _results(sandbox):
    return sandbox / "antarctica" / "results"


def _bank_a_series(sandbox, stem, years):
    r"""What a prior ISMIP7_OUTPUT=1 pass leaves beside the run's own files."""
    R = _results(sandbox)
    for y in years:
        (R / f"{stem}_ismip7_annual_{y}.h5").write_bytes(b"")
    (R / f"{stem}_ismip7_scalars.csv").write_text("year,lim\n")


def test_ocx_is_retried_from_its_own_saved_state(sandbox):
    r"""A wall stop at 2011 is resumed and reaches the 2026 target. Under the
    old `cold` classification the matrix stopped after one attempt and said
    the driver ignores ISMIP7_RESTART, which it no longer does."""
    out = _run(sandbox, FAKE_STEM="ocx_32000", FAKE_START_YEAR="1979",
               FAKE_STOP_YEAR="2011", FAKE_TARGET="2026")
    assert "COMPLETE at 2026" in out, out
    assert "not retried" not in out
    assert _attempts(sandbox) == ["cold", "restart"], out


def test_a_core_that_stops_advancing_is_still_left_alone(sandbox):
    r"""Retrying is gated on progress, not on the kind: a resumed attempt that
    reaches the same year as before gives up rather than looping."""
    out = _run(sandbox, FAKE_STEM="ocx_32000", FAKE_START_YEAR="1979",
               FAKE_STOP_YEAR="2011", FAKE_TARGET="2011")
    assert "no progress" in out, out
    assert _attempts(sandbox) == ["cold", "restart"], out


def test_a_fresh_rerun_archives_the_ismip7_series_too(sandbox):
    r"""FRESH=1 over a core that banked an ISMIP7 series must leave nothing
    behind that blocks the cold start. The series is moved, not deleted, so
    the operator keeps the audit trail."""
    stem = "ocx_32000"
    _bank_a_series(sandbox, stem, [2015, 2016, 2017])
    # a prior pass' own run products, which archive_stale already handled
    (_results(sandbox) / f"{stem}_timeseries.csv").write_text("year,vaf\n1979.0,0.0\n")
    (_results(sandbox) / f"{stem}_final.h5").write_bytes(b"")

    out = _run(sandbox, FAKE_STEM=stem, FAKE_START_YEAR="1979",
               FAKE_STOP_YEAR="2026", FAKE_TARGET="2026")

    R = _results(sandbox)
    assert list(R.glob(f"{stem}_ismip7_annual_*.h5")) == [], out
    assert not (R / f"{stem}_ismip7_scalars.csv").exists(), out
    archived = sorted(p.name for d in R.glob("archive_stale_*") for p in d.iterdir())
    assert f"{stem}_ismip7_scalars.csv" in archived, archived
    assert sum(1 for n in archived if "_ismip7_annual_" in n) == 3, archived
    assert "3 ISMIP7 year(s)" in out, out
    assert "COMPLETE at 2026" in out, out


def test_a_partially_complete_core_resumes_on_the_first_attempt(sandbox):
    r"""A core left short by a previous invocation keeps its output (nothing
    is stale), so attempt 1 must resume from its own checkpoint rather than
    cold-start, which would rewrite the years it banked."""
    import h5py

    stem = "ocx_32000"
    R = _results(sandbox)
    (R / f"{stem}_timeseries.csv").write_text("year,vaf\n1979.0,0.0\n2011.0,0.0\n")
    with h5py.File(R / f"{stem}_final.h5", "w") as h:
        h["/"].attrs["t_yr"] = 2011.0
        h["/"].attrs["stalled"] = 0
    _bank_a_series(sandbox, stem, [1979, 1980])

    out = _run(sandbox, FRESH="", FAKE_STEM=stem, FAKE_START_YEAR="1979",
               FAKE_STOP_YEAR="2011", FAKE_TARGET="2026")

    assert _attempts(sandbox) == ["restart"], out
    assert "resuming ocx_32000_final.h5 at t=2011" in out, out
    assert "COMPLETE at 2026" in out, out
    # nothing was archived: the output was current, so it is still in place
    assert list(R.glob("archive_stale_*")) == [], out
    assert len(list(R.glob(f"{stem}_ismip7_annual_*.h5"))) == 2, out


def test_a_core_with_no_state_of_its_own_still_cold_starts(sandbox):
    r"""The first ever run of a core has nothing to resume, and a branch core
    must reach its driver's own restart logic rather than be handed one."""
    out = _run(sandbox, FRESH="", FAKE_STEM="ocx_32000", FAKE_START_YEAR="1979",
               FAKE_STOP_YEAR="2026", FAKE_TARGET="2026")
    assert _attempts(sandbox) == ["cold"], out
    assert "starting (target 2026)" in out, out
