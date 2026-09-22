r"""The melt calibration follows the run's geometry space.

The forward melts DG0 cells with a cell slope and no cap; the earlier
calibrations fitted K on CG1 nodes with a capped slope, and the forward then
applied about 1.6 times the melt its K was fitted to (issue #30). The
calibration now reads ``ISMIP7_GEOMETRY_SPACE`` like the forward, records the
space it fitted on, and the loader says so once when a run melts on another.
"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("firedrake")
pytest.importorskip("rasterio")
pytest.importorskip("icepack")

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SCRIPTS = os.path.join(_REPO, "antarctica", "scripts")


@pytest.fixture
def calibrate_melt(monkeypatch):
    def _import(**env):
        for key in ("ISMIP7_GEOMETRY_SPACE", "ISMIP7_SIN_ALPHA_CAP",
                    "ISMIP7_MELT_OBS_CSV", "ISMIP7_K_OUT"):
            monkeypatch.delenv(key, raising=False)
        for key, value in env.items():
            monkeypatch.setenv(key, str(value))
        monkeypatch.syspath_prepend(_SCRIPTS)
        sys.modules.pop("calibrate_melt", None)
        import calibrate_melt
        return calibrate_melt
    yield _import
    sys.modules.pop("calibrate_melt", None)


def test_a_dg0_calibration_fits_the_slope_the_forward_melts_with(calibrate_melt):
    cm = calibrate_melt(ISMIP7_GEOMETRY_SPACE="dg0")
    assert cm.GEOMETRY == "dg0"
    assert cm.SIN_ALPHA_CAP == float("inf")


def test_a_cg1_calibration_keeps_its_cap(calibrate_melt):
    cm = calibrate_melt(ISMIP7_GEOMETRY_SPACE="cg1")
    assert cm.GEOMETRY == "cg1"
    assert cm.SIN_ALPHA_CAP == pytest.approx(5e-3)


def test_a_named_cap_applies_on_either_geometry(calibrate_melt):
    cm = calibrate_melt(ISMIP7_GEOMETRY_SPACE="dg0", ISMIP7_SIN_ALPHA_CAP="5e-3")
    assert cm.SIN_ALPHA_CAP == pytest.approx(5e-3)


def test_the_per_basin_fit_is_the_ratio_of_totals(calibrate_melt):
    cm = calibrate_melt()
    rho = float(cm._RHO_I)
    # Two basins of two dofs each and a third no dof melts in. Areas in m^2
    # and melt in m/yr chosen so the K = 1 totals are 1 and 2 Gt/yr.
    basin = np.array([1, 1, 2, 2, 5])
    area = np.array([1e9, 1e9, 1e9, 1e9, 1e9]) / rho * 1e3
    melt_1 = np.array([0.5, 0.5, 1.0, 1.0, 0.0])
    bids = np.array([1, 2, 3])
    M_obs = np.array([2.0, 1.0, 4.0])
    sigma = np.array([1.0, 2.0, 1.0])
    fit = cm.fit_per_basin_K(basin, melt_1, area, bids, M_obs, sigma)
    assert fit["M_1"] == pytest.approx([1.0, 2.0, 0.0])
    # K_b = M_obs / M_1 per basin, NaN where nothing melts
    assert fit["K_basin"][:2] == pytest.approx([2.0, 0.5])
    assert np.isnan(fit["K_basin"][2])
    # the Term-1 weighted scalar: sum(M_obs M_1 w) / sum(M_1^2 w), w = 1/sigma^2
    w = 1.0 / sigma ** 2
    M_1 = np.array([1.0, 2.0, 0.0])
    assert fit["K_star"] == pytest.approx(np.sum(M_obs * M_1 * w) / np.sum(M_1 ** 2 * w))
    # the total-match scalar: 7 / 3
    assert fit["K_total"] == pytest.approx(7.0 / 3.0)


@pytest.fixture
def basin_grid(tmp_path, monkeypatch):
    r"""A two-basin IMBIE-style grid under a data root, and a mesh of two
    points, one in each basin."""
    xr = pytest.importorskip("xarray")
    d = tmp_path / "parameterisations" / "ocean" / "imbie2"
    d.mkdir(parents=True)
    x = np.array([0.0, 8000.0])
    y = np.array([0.0, 8000.0])
    bn = np.array([[1, 1], [2, 2]], dtype="float32")
    xr.Dataset({"basinNumber": (("y", "x"), bn)}, coords={"x": x, "y": y}).to_netcdf(
        d / "basin_numbers_ismip8km_v2.nc")
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(tmp_path))
    return np.array([0.0, 0.0]), np.array([0.0, 8000.0])


def _npz(path, **extra):
    np.savez(path, basin_ids=np.array([1, 2]), K_basin=np.array([1e-4, 2e-4]), **extra)
    return str(path)


def test_the_loader_says_once_when_the_k_was_fitted_on_another_geometry(
        basin_grid, tmp_path, monkeypatch, capsys):
    import icepack2_tools.forcing as forcing
    monkeypatch.setattr(forcing, "_GEOMETRY_SPACE_WARNED", False)
    monkeypatch.setenv("ISMIP7_GEOMETRY_SPACE", "dg0")
    mx, my = basin_grid
    path = _npz(tmp_path / "K.npz", geometry_space="cg1")
    K = forcing.load_K_per_basin(path, mx, my)
    assert K == pytest.approx([1e-4, 2e-4])
    forcing.load_K_per_basin(path, mx, my)
    out = capsys.readouterr().out
    assert out.count("calibrated on cg1 geometry and this run melts on dg0") == 1


def test_a_file_without_the_tag_was_fitted_on_nodes(basin_grid, tmp_path, monkeypatch, capsys):
    import icepack2_tools.forcing as forcing
    monkeypatch.setattr(forcing, "_GEOMETRY_SPACE_WARNED", False)
    monkeypatch.setenv("ISMIP7_GEOMETRY_SPACE", "dg0")
    mx, my = basin_grid
    forcing.load_K_per_basin(_npz(tmp_path / "K.npz"), mx, my)
    assert "calibrated on cg1 geometry" in capsys.readouterr().out


def test_a_matching_geometry_is_silent(basin_grid, tmp_path, monkeypatch, capsys):
    import icepack2_tools.forcing as forcing
    monkeypatch.setattr(forcing, "_GEOMETRY_SPACE_WARNED", False)
    monkeypatch.setenv("ISMIP7_GEOMETRY_SPACE", "dg0")
    mx, my = basin_grid
    forcing.load_K_per_basin(_npz(tmp_path / "K.npz", geometry_space="dg0",
                                  sin_alpha_cap=float("inf")), mx, my)
    assert "WARNING" not in capsys.readouterr().out
