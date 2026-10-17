r"""The core report states a check it could not run as unjudged.

check_ismip6_track.py and compare_ismip6.py exit 2 when they cannot judge a
run, and compare_ismip6.py does so on every machine without the ISMIP6
ensemble archive (neither Quartz nor the Mac holds one). The report read
every non-zero status as a failure, so an unjudged run went on record as off
track or outside the envelope.
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

import core_report  # noqa: E402


def test_a_pass_and_a_failure_keep_their_words():
    assert core_report.verdict(0, "ON TRACK", "OFF TRACK") == "ON TRACK"
    assert core_report.verdict(1, "ON TRACK", "OFF TRACK") == "OFF TRACK"
    assert core_report.verdict(0, "inside envelope", "outside envelope") == "inside envelope"
    assert core_report.verdict(1, "inside envelope", "outside envelope") == "outside envelope"


def test_a_check_that_could_not_run_is_neither():
    words = core_report.verdict(2, "inside envelope", "outside envelope")
    assert "outside" not in words and "inside" not in words
    assert words.startswith("not judged") and "exit 2" in words


SCRIPTS = os.path.join(REPO, "antarctica", "scripts")
HEADER = ("year,mass_gt,vaf_mm_sle,smb_gtyr,melt_gtyr,outflux_gtyr,"
          "calv_gt,resid_gt\n")


def run(script, *args):
    import subprocess
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args],
                          capture_output=True, text=True)


def test_track_check_on_an_empty_timeseries_exits_2(tmp_path):
    csv_fn = tmp_path / "core01_ssp126_timeseries.csv"
    csv_fn.write_text(HEADER)
    r = run("check_ismip6_track.py", str(csv_fn))
    assert r.returncode == 2
    assert "empty timeseries" in r.stderr


def test_track_check_on_a_missing_timeseries_exits_2(tmp_path):
    r = run("check_ismip6_track.py", str(tmp_path / "absent_timeseries.csv"))
    assert r.returncode == 2
    assert "absent_timeseries.csv" in r.stderr


def test_ensemble_check_on_a_missing_control_exits_2(tmp_path):
    proj = tmp_path / "core01_ssp126_timeseries.csv"
    proj.write_text(HEADER + "2015.0,1,58000,2400,1100,1200,0,0\n"
                             "2016.0,1,57999,2400,1100,1200,0,0\n")
    r = run("compare_ismip6.py", str(proj), str(tmp_path / "ctrl_timeseries.csv"),
            "--ismip6-dir", str(tmp_path))
    assert r.returncode == 2
    assert "ctrl_timeseries.csv" in r.stdout


def test_ensemble_check_on_an_empty_projection_exits_2(tmp_path):
    proj = tmp_path / "core01_ssp126_timeseries.csv"
    proj.write_text(HEADER)
    r = run("compare_ismip6.py", str(proj), str(proj), "--ismip6-dir", str(tmp_path))
    assert r.returncode == 2
