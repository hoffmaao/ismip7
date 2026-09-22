r"""The protocol's per-basin adjustment is a thermal-forcing offset at one
toolbox K, fitted on the forward's melt path and stamped onto any mesh."""
import importlib.util
import os

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _calibrate_deltaT():
    spec = importlib.util.spec_from_file_location(
        "calibrate_deltaT",
        os.path.join(REPO, "antarctica", "scripts", "calibrate_deltaT.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _synthetic(K=8.5e-5, true_dT=(0.3, -0.5)):
    from icepack2_tools.forcing import quadratic_mixed_slope, _RHO_I
    rng = np.random.default_rng(0)
    n = 400
    tf = rng.uniform(0.5, 2.5, n)
    sal = np.full(n, 34.5)
    sin_a = np.full(n, 5.115e-3)
    area = np.full(n, 4.0e6)
    basin = np.repeat([3, 9], n // 2)
    floating = np.ones(n, bool)
    bids = np.array([3, 9])
    M_obs = np.zeros(2)
    for i, (bid, d) in enumerate(zip(bids, true_dT)):
        sel = basin == bid
        m = quadratic_mixed_slope(tf[sel] + d, sal[sel], sin_a[sel], K=K)
        M_obs[i] = float((m * area[sel]).sum()) * float(_RHO_I) / 1e12
    return dict(tf=tf, sal=sal, sin_a=sin_a, K=K, floating=floating,
                area=area, basin=basin, bids=bids, M_obs=M_obs)


def test_the_fit_recovers_a_known_offset():
    cd = _calibrate_deltaT()
    d = _synthetic()
    dT, resid, sens, flagged = cd.fit_deltaT(
        d["tf"], d["sal"], d["sin_a"], d["K"], d["floating"], d["area"],
        d["basin"], d["bids"], d["M_obs"])
    assert np.allclose(dT, (0.3, -0.5), atol=2e-4)
    assert np.all(np.abs(resid) < 1e-3)  # Gt/yr, the root tolerance in dT
    assert np.all(sens > 0)
    assert flagged == []


def test_an_unreachable_basin_takes_the_window_end_and_is_flagged():
    cd = _calibrate_deltaT()
    d = _synthetic()
    d["M_obs"][0] *= 50.0
    dT, resid, sens, flagged = cd.fit_deltaT(
        d["tf"], d["sal"], d["sin_a"], d["K"], d["floating"], d["area"],
        d["basin"], d["bids"], d["M_obs"])
    assert dT[0] == cd.DT_WINDOW[1]
    assert [b for b, _ in flagged] == [3]


def test_the_offset_is_stamped_by_basin(tmp_path, monkeypatch):
    xr = pytest.importorskip("xarray")
    from icepack2_tools.forcing import load_deltaT_per_basin
    # a 4 x 4 basin grid: basin 3 on the left half, basin 9 on the right
    x = np.array([0.0, 1.0, 2.0, 3.0])
    y = np.array([0.0, 1.0, 2.0, 3.0])
    bn = np.where(np.arange(4)[None, :] < 2, 3, 9).repeat(4, axis=0)
    imbie = tmp_path / "basins.nc"
    xr.Dataset({"basinNumber": (("y", "x"), bn)},
               coords={"x": x, "y": y}).to_netcdf(imbie)
    npz = tmp_path / "deltaT.npz"
    np.savez(npz, basin_ids=np.array([3, 9]), deltaT_basin=np.array([0.3, -0.5]),
             K=8.5e-5, melt_slope="ant", imbie2_nc=str(imbie))
    monkeypatch.setenv("ISMIP7_MELT_SLOPE", "ant")
    mx = np.array([0.4, 2.6, 0.2, 2.9])
    my = np.array([0.1, 0.1, 2.9, 2.9])
    field, K = load_deltaT_per_basin(str(npz), mx, my)
    assert K == 8.5e-5
    assert np.allclose(field, [0.3, -0.5, 0.3, -0.5])
