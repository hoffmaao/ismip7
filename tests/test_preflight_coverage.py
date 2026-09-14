r"""What the preflight gate calls a missing input.

The reader bridges exactly one year past the end of a series (CESM2-WACCM's
atmosphere stops at 2299 and the empty 2300 files were removed, discussion
#8), so a tree whose LAST year is absent runs correctly. The gate has to agree:
calling that a missing input reports cores 5 and 7 as BLOCKED for the
production tree they are meant to run against. Any other hole stays an error,
because the reader raises on it.

``atm_years`` reads only filenames, so the synthetic trees here are empty
files in the layout ``atmosphere_path`` resolves. The shared mesh/MAP checks
are stubbed out so the status reflects the forcing coverage alone; everything
else, including the printed status line, is the shipped code.
"""
import importlib
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

pytest.importorskip("firedrake")

ESM, SCENARIO = "CESM2-WACCM", "ssp585"
CORE_7 = "core  7"


def _tree(root, years, esm=ESM, scenario=SCENARIO):
    r"""Empty acabf-anomaly files for `years`, in the layout the reader uses."""
    for var in ("acabf-anomaly", "acabf"):
        d = os.path.join(root, esm, scenario, "SDBN1-8000m", var, "v2")
        os.makedirs(d, exist_ok=True)
        for y in years:
            head = f"{var}_AIS_{esm}_{scenario}_SDBN1-8000m_v2_{y}.nc"
            open(os.path.join(d, head), "wb").close()


def _run(monkeypatch, tmp_path, years):
    r"""preflight.main() over a tree covering `years`, returning its output."""
    root = str(tmp_path / "ISMIP7" / "AIS")
    os.makedirs(root, exist_ok=True)
    _tree(root, years)
    monkeypatch.setenv("ISMIP7_DATA_ROOT", root)
    sys.modules.pop("preflight", None)
    preflight = importlib.import_module("preflight")
    # isolate the forcing-coverage verdict from the mesh/MAP/ocean checks
    monkeypatch.setattr(preflight, "shared_missing", lambda warn: [])
    monkeypatch.setattr(preflight, "racmo_ok", lambda: True)
    monkeypatch.setattr(preflight, "oi_ok", lambda: True)
    monkeypatch.setattr(preflight, "ocean_cover", lambda e, s: (2015, 2300))
    monkeypatch.setattr(preflight, "pool_status",
                        lambda *a, **k: ("ok", ""))
    preflight.main()
    return None


def _core_line(capsys, core=CORE_7):
    line = [ln for ln in capsys.readouterr().out.splitlines() if core in ln]
    assert len(line) == 1, line
    return line[0]


def test_an_absent_final_year_is_a_note_not_a_blocker(monkeypatch, tmp_path, capsys):
    r"""CESM2-WACCM ssp585 covers 2015-2300 and the tree stops at 2299: the
    reader persists 2299 for 2300, so the core is READY."""
    _run(monkeypatch, tmp_path, range(2015, 2300))       # 2015..2299
    line = _core_line(capsys)
    assert "READY" in line, line
    assert "BLOCKED" not in line
    assert "2300 is absent" in line and "persists 2299" in line


def test_a_hole_inside_the_series_still_blocks(monkeypatch, tmp_path, capsys):
    r"""The reader raises on an interior gap, so the gate must not clear it."""
    years = [y for y in range(2015, 2301) if y != 2100]
    _run(monkeypatch, tmp_path, years)
    line = _core_line(capsys)
    assert "BLOCKED" in line, line
    assert "2100..2100" in line


def test_a_short_tree_still_blocks(monkeypatch, tmp_path, capsys):
    r"""Two years short is not the one-year bridge; it is a short download."""
    _run(monkeypatch, tmp_path, range(2015, 2299))       # 2015..2298
    line = _core_line(capsys)
    assert "BLOCKED" in line, line


def test_a_complete_tree_is_ready_without_a_note(monkeypatch, tmp_path, capsys):
    _run(monkeypatch, tmp_path, range(2015, 2301))       # 2015..2300
    line = _core_line(capsys)
    assert "READY" in line, line
    assert "absent" not in line
