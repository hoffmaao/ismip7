r"""A K file records the draft slope it was fitted against, and the forward
says so when that disagrees with what it applies.

Melt is linear in ``sin(alpha)``, so a per-basin K is only valid for the slope
field it was calibrated on. ``calibrate_melt.py`` caps the slope at
``ISMIP7_SIN_ALPHA_CAP`` and stamps that value into the npz;
``forcing.compute_sin_alpha``, which the forward calls every step, applies no
cap. At the reference state on the Úa 2 km mesh the forward's own uncapped
operator integrates about 4293 Gt/yr against the 1067.4 Gt/yr the K was fitted
to, and capping that operator lands near 1028. The mismatch used to be
invisible: both halves ran, and nothing said they disagreed.

These tests pin the warning, and the physics stays as it is. Which side should
move is a science decision tied to the unsettled upstream local-slope
question, so the forward's numbers are deliberately unchanged.

Every test goes through ``forcing.load_K_per_basin``, the loader the forward
calls, over a synthetic IMBIE2 basin file in the layout it resolves under
``ISMIP7_DATA_ROOT``. Serial, one small npz, no mesh.
"""

import os

import numpy as np
import pytest

pytest.importorskip("firedrake")
xr = pytest.importorskip("xarray")

from icepack2_tools import forcing                                # noqa: E402

CAP = 5e-3


@pytest.fixture(autouse=True)
def _rearm():
    r"""The warning fires once per process, so re-arm it around each test."""
    before = forcing._SLOPE_CAP_WARNED
    forcing._SLOPE_CAP_WARNED = False
    yield
    forcing._SLOPE_CAP_WARNED = before


@pytest.fixture
def basins(tmp_path, monkeypatch):
    r"""A 2 by 2 IMBIE2 grid with basins 1 and 2 side by side."""
    root = tmp_path / "ISMIP7" / "AIS"
    d = root / "parameterisations" / "ocean" / "imbie2"
    os.makedirs(d)
    xr.Dataset(
        {"basinNumber": (("y", "x"), np.array([[1, 2], [1, 2]]))},
        coords={"x": [0.0, 8000.0], "y": [0.0, 8000.0]},
    ).to_netcdf(d / "basin_numbers_ismip8km_v2.nc")
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(root))


def _npz(path, **extra):
    np.savez(path, basin_ids=np.array([1, 2]),
             K_basin=np.array([1e-4, 2e-4]), **extra)
    return str(path)


def _load(path):
    return forcing.load_K_per_basin(path, [0.0, 8000.0], [0.0, 0.0])


def test_a_capped_calibration_is_announced(basins, tmp_path, capsys):
    path = _npz(tmp_path / "calibrated_K_per_basin_2000.npz",
                sin_alpha_cap=CAP)
    K = _load(path)
    said = capsys.readouterr().out

    np.testing.assert_allclose(K, [1e-4, 2e-4])
    assert "calibrated_K_per_basin_2000.npz" in said
    assert "0.005" in said
    # The reader has to be able to find out what to do about it.
    assert "FORWARD_RUN_READINESS" in said


def test_it_says_so_once(basins, tmp_path, capsys):
    path = _npz(tmp_path / "a.npz", sin_alpha_cap=CAP)
    _load(path)
    assert capsys.readouterr().out.strip()

    # A chained resume or a driver that reloads must not turn this into a
    # per-step banner.
    _load(path)
    assert capsys.readouterr().out == ""


def test_a_file_without_the_stamp_is_not_flagged(basins, tmp_path, capsys):
    r"""An older npz carries no cap, so there is nothing to compare and the
    loader stays quiet and still returns the K field."""
    path = _npz(tmp_path / "old.npz")
    K = _load(path)

    np.testing.assert_allclose(K, [1e-4, 2e-4])
    assert capsys.readouterr().out == ""
