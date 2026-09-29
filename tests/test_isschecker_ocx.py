r"""``isschecker_ocx.py``: isschecker 0.5.1 with core 11 checked (issue #18).

The patch adds the ``ocx`` experiment and accepts ``ERA5`` in field 5 of ocx
files, and leaves every other file to the release's own checks. Most of it
is exercised on a stand-in for the checker module, so it runs without
isschecker; the last test holds the installed release to the names the patch
reaches into.
"""
import os
import sys
import types

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "antarctica", "scripts"))
import isschecker_ocx  # noqa: E402

OCX_FILE = "lithk_AIS_RICE_icepack2_m001_ERA5_f001_ocx_C011_2003-2025.nc"
CTRL_FILE = "lithk_AIS_RICE_icepack2_m001_ERA5_f001_ctrl_C009_2015-2300.nc"


def _checker(version="0.5.1", experiments=("historical", "ctrl"), names=("CESM2-WACCM",)):
    r"""A stand-in for ``isschecker.checker``, recording the field 5 names
    each naming check saw."""
    seen = []

    def load(*args, **kwargs):
        return [{"experiment": e, "start_year_min": 2015, "start_year_max": 2015,
                 "end_year": 2300, "duration": 286} for e in experiments]

    checker = types.SimpleNamespace(
        __version__=version, VALID_ESM_NAMES=set(names),
        ISMIP7_FILENAME_PARTS=10, ISMIP7_FILENAME_EXPERIMENT_IDX=7,
        _load_experiments_csv=load, _describe_version=lambda: version)

    def check_naming(reporter, file_name, *args):
        seen.append(set(checker.VALID_ESM_NAMES))
        return file_name

    checker._check_naming = check_naming
    return checker, seen


def test_the_ocx_experiment_is_the_protocol_overview_row():
    checker, _ = _checker()
    line = isschecker_ocx.patch(checker)
    rows = checker._load_experiments_csv()
    assert [r["experiment"] for r in rows] == ["historical", "ctrl", "ocx"]
    assert rows[-1] == {"experiment": "ocx", "start_year_min": 1990, "start_year_max": 2015,
                        "end_year": 2025, "duration": -1}
    assert line == "0.5.1 " + isschecker_ocx.PATCH_NOTE


def test_era_is_accepted_in_field_5_of_ocx_files_only():
    checker, seen = _checker()
    isschecker_ocx.patch(checker)
    assert checker._check_naming(None, OCX_FILE, "AIS") == OCX_FILE
    checker._check_naming(None, CTRL_FILE, "AIS")
    assert "ERA5" in seen[0] and "ERA5" not in seen[1]
    assert checker.VALID_ESM_NAMES == {"CESM2-WACCM"}


def test_a_name_the_release_already_accepts_stays_accepted():
    checker, _ = _checker(names=("CESM2-WACCM", "ERA5"))
    isschecker_ocx.patch(checker)
    checker._check_naming(None, OCX_FILE)
    assert "ERA5" in checker.VALID_ESM_NAMES


def test_another_release_is_refused():
    checker, _ = _checker(version="0.5.2")
    with pytest.raises(SystemExit, match="written against 0.5.1"):
        isschecker_ocx.patch(checker)


def test_a_release_with_its_own_ocx_experiment_retires_the_patch():
    checker, _ = _checker(experiments=("historical", "ocx"))
    isschecker_ocx.patch(checker)
    with pytest.raises(SystemExit, match="ocx experiment of its own"):
        checker._load_experiments_csv()


def test_the_installed_release_takes_the_patch(monkeypatch):
    checker = pytest.importorskip("isschecker.checker")
    if checker.__version__ != isschecker_ocx.PATCHED_RELEASE:
        pytest.skip(f"isschecker {checker.__version__} installed; the script itself refuses it")
    for name in ("_load_experiments_csv", "_check_naming"):
        monkeypatch.setattr(checker, name, getattr(checker, name))
    for name in ("_describe_version", "_parse_args", "run_checker"):
        assert callable(getattr(checker, name))
    isschecker_ocx.patch(checker)
    assert checker._load_experiments_csv()[-1] == isschecker_ocx.OCX_EXPERIMENT
    assert "ERA5" not in checker.VALID_ESM_NAMES
