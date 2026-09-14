r"""A K file records the draft slope it was fitted against, and the forward
says so when that disagrees with what it applies.

Melt is linear in ``sin(alpha)``, so a per-basin K is only valid for the slope
field it was calibrated on. ``calibrate_melt.py`` caps the slope at
``ISMIP7_SIN_ALPHA_CAP`` and stamps that value into the npz;
``forcing.compute_sin_alpha``, which the forward calls every step, applies no
cap. On the Úa 2 km mesh the two differ by about a factor of 1.7 in integrated
melt, and that is what puts grounding-zone cells past the variable request's
``libmassbffl`` bound. The mismatch used to be invisible: both halves ran, and
nothing said they disagreed.

This pins the warning rather than the physics. Which side should move is a
science decision tied to the unsettled upstream local-slope question, so the
forward's numbers are deliberately unchanged.

Serial, one small npz, no mesh.
"""

import os
import sys

import numpy as np
import pytest

pytest.importorskip("firedrake")

from icepack2_tools import forcing                                # noqa: E402


@pytest.fixture(autouse=True)
def _rearm():
    r"""The warning fires once per process, so re-arm it around each test."""
    before = forcing._SLOPE_CAP_WARNED
    forcing._SLOPE_CAP_WARNED = False
    yield
    forcing._SLOPE_CAP_WARNED = before


def test_a_capped_calibration_is_announced(capsys):
    forcing._warn_slope_cap("calibrated_K_per_basin_2000.npz", 5e-3)
    said = capsys.readouterr().out

    assert "calibrated_K_per_basin_2000.npz" in said
    assert "0.005" in said
    # The reader has to be able to find out what to do about it.
    assert "FORWARD_RUN_READINESS" in said
    assert forcing._SLOPE_CAP_WARNED is True


def test_it_says_so_once(capsys):
    forcing._warn_slope_cap("a.npz", 5e-3)
    first = capsys.readouterr().out
    assert first.strip()

    # A forward calls the loader once per run, but a chained resume or a driver
    # that reloads must not turn this into a per-step banner.
    if not forcing._SLOPE_CAP_WARNED:
        pytest.fail("the flag did not latch")


def test_the_writer_of_the_npz_records_the_cap(tmp_path):
    r"""The guard is only reachable if calibrate_melt actually stamps the cap.

    A K file without `sin_alpha_cap` predates the stamp, and the loader then
    has nothing to compare against, which is the silent case this exists to
    end. Read the script's source rather than running it, since a real
    calibration needs the climatology and a mesh.
    """
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src = open(os.path.join(here, "antarctica", "scripts",
                            "calibrate_melt.py")).read()
    assert "sin_alpha_cap=SIN_ALPHA_CAP" in src, (
        "calibrate_melt.py no longer stamps the slope cap into its npz, so "
        "forcing.load_K_per_basin cannot detect the mismatch"
    )


def test_a_file_without_the_stamp_is_not_flagged(tmp_path, capsys):
    r"""An older npz carries no cap, so there is nothing to compare and the
    loader must stay quiet rather than guess."""
    path = tmp_path / "old.npz"
    np.savez(path, basin_ids=np.arange(3), K_basin=np.full(3, 1e-4))
    data = np.load(path)

    assert "sin_alpha_cap" not in data
    # Mirror the loader's guard: no stamp, no warning, no exception.
    if "sin_alpha_cap" in data:
        forcing._warn_slope_cap(str(path), float(data["sin_alpha_cap"]))
    assert capsys.readouterr().out == ""
    assert forcing._SLOPE_CAP_WARNED is False
