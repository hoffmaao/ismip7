r"""The melt slope follows ISMIP7_MELT_SLOPE.

The ISMIP7 reference example is the quadratic law with one constant mean
Antarctic slope; the local slope is its other option. The forward's slope,
the calibration and the K loader read one knob, and a K fitted under one
convention announces itself once under the other.
"""
import os

import numpy as np
import pytest

fd = pytest.importorskip("firedrake")
xr = pytest.importorskip("xarray")

import icepack2_tools.forcing as forcing  # noqa: E402


@pytest.fixture(autouse=True)
def _clean(monkeypatch):
    for k in ("ISMIP7_MELT_SLOPE", "ISMIP7_SIN_ALPHA_ANT", "ISMIP7_GEOMETRY_SPACE"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setattr(forcing, "_MELT_SLOPE_WARNED", False)
    monkeypatch.setattr(forcing, "_SLOPE_CAP_WARNED", False)
    monkeypatch.setattr(forcing, "_GEOMETRY_SPACE_WARNED", False)


def _ctx(degree):
    mesh = fd.RectangleMesh(4, 4, 8000.0, 8000.0)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    Q_g = fd.FunctionSpace(mesh, "DG", 0) if degree == 0 else Q
    x, y = fd.SpatialCoordinate(mesh)
    h = fd.Function(Q_g).interpolate(500.0 + 0.05 * x)      # a sloping draft
    s = fd.Function(Q_g).interpolate(50.0 + 0.005 * x)
    return {"Q": Q, "V": V, "Q_g": Q_g, "h": h, "s": s}


@pytest.mark.parametrize("degree", [0, 1])
def test_the_default_is_one_constant_slope_on_every_dof(degree):
    ctx = _ctx(degree)
    sin_a = forcing.compute_sin_alpha(ctx)
    assert sin_a.shape == ctx["h"].dat.data_ro.shape
    assert np.allclose(sin_a, forcing.SIN_ALPHA_ANT_DEFAULT)


def test_the_constant_is_a_knob(monkeypatch):
    monkeypatch.setenv("ISMIP7_SIN_ALPHA_ANT", "4e-3")
    assert np.allclose(forcing.compute_sin_alpha(_ctx(0)), 4e-3)


@pytest.mark.parametrize("degree", [0, 1])
def test_local_is_the_slope_of_this_draft(monkeypatch, degree):
    monkeypatch.setenv("ISMIP7_MELT_SLOPE", "local")
    sin_a = forcing.compute_sin_alpha(_ctx(degree))
    # draft = s - h has slope -0.045 in x, so sin(alpha) = 0.045 / sqrt(1 + 0.045^2).
    # On CG1 the projected gradient is exact; the DG0 path lifts cell means to
    # CG1 first, which is exact inside and biased along the boundary.
    expect = 0.045 / np.sqrt(1.0 + 0.045 ** 2)
    if degree == 1:
        assert np.allclose(sin_a, expect, rtol=1e-6)
    else:
        assert np.median(sin_a) == pytest.approx(expect, rel=0.15)
        assert not np.allclose(sin_a, sin_a[0])


def test_an_unknown_convention_is_refused(monkeypatch):
    monkeypatch.setenv("ISMIP7_MELT_SLOPE", "mean")
    with pytest.raises(ValueError):
        forcing.melt_slope()


@pytest.fixture
def basins(tmp_path, monkeypatch):
    root = tmp_path / "ISMIP7" / "AIS"
    d = root / "parameterisations" / "ocean" / "imbie2"
    os.makedirs(d)
    xr.Dataset({"basinNumber": (("y", "x"), np.array([[1, 2], [1, 2]]))},
               coords={"x": [0.0, 8000.0], "y": [0.0, 8000.0]}).to_netcdf(
        d / "basin_numbers_ismip8km_v2.nc")
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(root))


def _npz(path, **extra):
    np.savez(path, basin_ids=np.array([1, 2]), K_basin=np.array([1e-4, 2e-4]), **extra)
    return str(path)


def _load(path):
    return forcing.load_K_per_basin(path, [0.0, 8000.0], [0.0, 0.0])


def test_a_k_fitted_under_the_other_convention_says_so_once(basins, tmp_path, capsys):
    path = _npz(tmp_path / "K.npz", melt_slope="local", sin_alpha_cap=5e-3)
    _load(path); _load(path)
    out = capsys.readouterr().out
    assert out.count("calibrated under ISMIP7_MELT_SLOPE=local and this run melts under ant") == 1
    # the cap warning belongs to the local convention and is not raised on top
    assert "capped at sin(alpha)" not in out


def test_a_file_without_the_entry_was_fitted_on_the_local_slope(basins, tmp_path, capsys):
    _load(_npz(tmp_path / "K.npz"))
    assert "ISMIP7_MELT_SLOPE=local" in capsys.readouterr().out


def test_a_matching_convention_is_silent(basins, tmp_path, capsys):
    K = _load(_npz(tmp_path / "K.npz", melt_slope="ant", sin_alpha_ant=2.9e-3,
                   sin_alpha_cap=float("inf"), geometry_space="dg0"))
    assert np.allclose(K, [1e-4, 2e-4])
    assert "WARNING" not in capsys.readouterr().out


def test_local_keeps_the_cap_warning(basins, tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("ISMIP7_MELT_SLOPE", "local")
    _load(_npz(tmp_path / "K.npz", melt_slope="local", sin_alpha_cap=5e-3))
    out = capsys.readouterr().out
    assert "capped at sin(alpha) = 0.005" in out and "ISMIP7_MELT_SLOPE=" not in out
