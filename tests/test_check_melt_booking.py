r"""check_melt_booking.py's arithmetic on a four-cell year (issue #136).

Cell 0 is a floating shelf that loses 1 m of ice across the mesh's exterior
boundary, cell 1 a marine cell holding no ice at either end that the inflow,
the frozen reference and snowfall feed, cell 2 an ice-free land cell the
reference feeds, cell 3 grounded ice on land that supplies the inflow. Pixel
0 covers cell 0, pixel 1 cell 1, pixel 2 cells 2 and 3, a thousand square
kilometres each. The same year is booked twice: as the forward did before issue #136,
the reference and the snowfall melted on cell 1 with only the inflow's share
split into lifmassbf, the land film melted with ligroundf booked into it, and
the outflux in no field; and as it does now.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "antarctica", "scripts"))
sparse = pytest.importorskip("scipy.sparse")
pytest.importorskip("firedrake")
cmb = pytest.importorskip("check_melt_booking")

from icepack2_tools.ismip7_output import net_reference, net_snowfall, split_melt  # noqa: E402

A = 1e9                               # m2: the budget gate works in Gt/yr
GT = A * 917.0 / 1e12                 # one metre of ice over one cell, in Gt
W = sparse.csr_matrix(np.array([[A, 0, 0, 0], [0, A, 0, 0], [0, 0, A, A]], dtype=float))


def _state(ice, floating):
    topg = np.array([-500.0, -500.0, 100.0, 100.0])
    lithk = np.where(ice, np.array([300.0, 0.0, 0.0, 200.0]), 0.0)
    orog = np.array([31.0, 0.0, 100.0, 300.0])
    ice = np.asarray(ice, dtype=float)
    fl = np.asarray(floating, dtype=float)
    return {"sftgif": ice, "sftflf": fl, "sftgrf": ice - fl, "topg": topg,
            "lithk": lithk, "orog": orog}


def _year(booking):
    r"""The start state and the year's end state and fields, in m/yr. The
    inflow is 2 m into cell 0 and 1 m into cell 1, from cell 3; now the
    reference on the land cell flows on to cell 3 instead of melting."""
    prev = _state([1, 0, 0, 1], [1, 0, 0, 0])
    cur = _state([1, 0, 0, 1], [1, 0, 0, 0])
    melt = np.array([-2.0, -6.0, -0.5, 0.0])     # the melt before the split
    smb = np.array([0.0, 1.0, 0.0, 0.0])
    ref = np.array([0.0, 4.0, 0.5, 0.0])
    dh = np.array([-1.0, 0.0, 0.0, -3.0])
    calv = np.zeros(4)
    lig = np.array([0.0, 0.0, 2.0, 0.0])
    none = np.array([False, True, True, False])
    if booking == "now":
        melt[2] = 0.0                            # the film takes no melt
        dh[3] = -2.5                             # its reference flows on instead
        lig[2] = 0.0                             # land is grounded to the booking
        calv[0] = -1.0                           # the outflux is calving
        melt, smb, ref = net_reference(melt, smb, ref, none)
        melt, smb = net_snowfall(melt, smb, none)
        lib, lif = split_melt(melt, cur["sftflf"] > 0.5)
    else:
        # the inflow's share of cell 1's melt, 1 m of its supply of 6
        lif = np.array([0.0, -1.0, 0.0, 0.0])
        lib = melt - lif
    cur.update(libmassbffl=lib, lifmassbf=lif, acabf=smb, acabf_correction=ref,
               dlithkdt=dh, licalvf=calv, ligroundf=lig)
    return prev, cur


def test_the_old_booking_shows_what_the_fields_missed():
    prev, cur = _year("before")
    r = cmb.year_booking(prev, cur, W)
    assert r["lifmassbf"] == pytest.approx(1.0 * GT)
    assert r["reference in sinks"] == pytest.approx(4.5 * GT)
    assert r["snowfall in sinks"] == pytest.approx(1.0 * GT)
    assert r["land melt"] == pytest.approx(0.5 * GT)
    assert r["land ligroundf"] == pytest.approx(2.0 * GT)
    assert r["left out"] == pytest.approx(5.5 * GT)        # cell 1's 5 m and the film
    assert r["reference on melting cells"] == pytest.approx(4.5 * GT)
    assert r["budget residual"] == pytest.approx(-1.0 * GT)   # the outflux
    assert r["split"] == pytest.approx(5.0)
    bad = cmb.failures([{"year": 2300, **r}], tol=1e-6)
    assert sorted(k for _, k, _ in bad) == ["budget residual", "land ligroundf", "land melt",
                                            "left out", "reference in sinks",
                                            "snowfall in sinks", "split"]


def test_the_current_booking_passes_the_gate():
    prev, cur = _year("now")
    r = cmb.year_booking(prev, cur, W)
    assert r["lifmassbf"] == pytest.approx(1.0 * GT)
    assert r["lifmassbf: no ice at either end"] == pytest.approx(1.0 * GT)
    assert r["lifmassbf: afloat at start"] == 0.0 and r["lifmassbf: the rest"] == 0.0
    assert r["reference in sinks"] == pytest.approx(0.0, abs=1e-15)
    assert r["snowfall in sinks"] == pytest.approx(0.0, abs=1e-15)
    assert r["land melt"] == 0.0 and r["land ligroundf"] == 0.0
    assert r["left out"] == 0.0 and r["refreezing left out"] == 0.0
    assert r["licalvf"] == pytest.approx(1.0 * GT)
    assert r["budget residual"] == pytest.approx(0.0, abs=1e-12)
    assert cmb.failures([{"year": 2300, **r}], tol=1e-6) == []


def test_a_forward_split_that_does_not_reproduce_fails():
    r"""Half of cell 1's front melt written back into libmassbffl: the melt
    is the same, and the file's lifmassbf is not what the rule gives, so the
    fill leaves half of it out."""
    prev, cur = _year("now")
    half = cur["lifmassbf"][1] / 2
    cur["lifmassbf"][1] -= half
    cur["libmassbffl"][1] += half
    r = cmb.year_booking(prev, cur, W)
    assert r["split"] == pytest.approx(0.5)
    assert r["left out"] == pytest.approx(0.5 * GT)
    assert sorted(k for _, k, _ in cmb.failures([{"year": 2300, **r}], tol=1e-6)) == [
        "left out", "split"]


def test_licalvf_is_held_to_the_forward_s_calving_and_outflux(tmp_path):
    r"""The timeseries books calving per step in Gt and the outflux as a
    rate: a year's sum of both is what licalvf carries, and a row repeated by
    a resumed link counts once. A year missing from the timeseries fails."""
    path = tmp_path / "ts.csv"
    rows = [(2300.25, 0.1, GT), (2300.5, 0.0, GT), (2300.5, 0.0, GT), (2300.75, 0.2, 0.0),
            (2301.0, 0.0, 2 * GT), (2301.25, 5.0, 0.0)]
    with open(path, "w") as f:
        f.write("year,outflux_gtyr,calv_gt\n")
        for t, calv, out in rows:
            f.write(f"{t:.6f},{out:.10f},{calv:.10f}\n")
    losses = cmb.timeseries_losses(str(path))
    assert losses[2300] == pytest.approx(0.3 + GT)          # 0.25 yr at GT, twice at 2 GT
    assert losses[2301] == pytest.approx(5.0)
    prev, cur = _year("now")
    r = {"year": 2300, **cmb.year_booking(prev, cur, W)}
    r["licalvf against the timeseries"] = r["licalvf"] - (losses[2300] - 0.3)
    assert cmb.failures([r], tol=1e-6) == []
    r["licalvf against the timeseries"] = r["licalvf"] - losses[2300]
    assert [k for _, k, _ in cmb.failures([r], tol=1e-6)] == ["licalvf against the timeseries"]
    r["licalvf against the timeseries"] = float("nan")
    assert [k for _, k, _ in cmb.failures([r], tol=1e-6)] == ["licalvf against the timeseries"]


def test_a_missing_series_is_an_input_error(tmp_path):
    assert cmb.main([str(tmp_path / "nothing_ismip7_annual.h5")]) == 2
