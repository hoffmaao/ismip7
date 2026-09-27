r"""check_melt_booking.py's arithmetic on a four-cell year (issue #136).

Cell 0 is a floating shelf, cell 1 a marine cell holding no ice at either end
that the inflow and the frozen reference feed, cell 2 an ice-free land cell
carrying a film, cell 3 grounded ice on land. Pixel 0 covers cell 0, pixel 1
cell 1, pixel 2 cells 2 and 3, one square kilometre each. The same year is
booked twice: as the forward did before issue #136, the reference melted on
cell 1 and the film on cell 2 melted with ligroundf booked into it, and as it
does now.
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

from icepack2_tools.ismip7_output import net_reference, split_front_melt  # noqa: E402

A = 1e6
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
    r"""The end state and the year's fields in m/yr."""
    prev = _state([1, 0, 0, 1], [1, 0, 0, 0])
    cur = _state([1, 0, 0, 1], [1, 0, 0, 0])
    melt = np.array([-2.0, -5.0, -0.5, 0.0])     # the melt before the split
    smb = np.zeros(4)
    ref = np.array([0.0, 4.0, 0.5, 0.0])
    dh = np.zeros(4)
    lig = np.array([0.0, 0.0, 2.0, 0.0])
    none = np.array([False, True, True, False])
    marine_none = np.array([False, True, False, False])
    if booking == "now":
        melt[2] = 0.0                            # the film takes no melt
        lig[2] = 0.0                             # land is grounded to the booking
        melt, smb, ref = net_reference(melt, smb, ref, none)
    lib, lif = split_front_melt(melt, smb, ref, dh, marine_none)
    cur.update(libmassbffl=lib, lifmassbf=lif, acabf=smb, acabf_correction=ref,
               dlithkdt=dh, ligroundf=lig)
    return prev, cur


def test_the_old_booking_shows_the_reference_and_the_land_film():
    prev, cur = _year("before")
    r = cmb.year_booking(prev, cur, W)
    assert r["lifmassbf"] == pytest.approx(1.0 * GT)       # 1 m of inflow of 5 supplied
    assert r["reference in sinks"] == pytest.approx(4.5 * GT)
    assert r["land melt"] == pytest.approx(0.5 * GT)
    assert r["land ligroundf"] == pytest.approx(2.0 * GT)
    assert r["left out"] == pytest.approx(4.5 * GT)        # cell 1's 4 m and the film
    assert r["left out: no ice at either end"] == pytest.approx(4.5 * GT)
    assert r["left out: afloat at start"] == pytest.approx(0.0)
    assert r["split"] == pytest.approx(0.0, abs=1e-15)
    bad = cmb.failures([{"year": 2300, **r}], tol=1e-6)
    assert sorted(k for _, k, _ in bad) == ["land ligroundf", "land melt", "reference in sinks"]


def test_the_current_booking_passes_the_gate():
    prev, cur = _year("now")
    r = cmb.year_booking(prev, cur, W)
    assert r["lifmassbf"] == pytest.approx(1.0 * GT)
    assert r["reference in sinks"] == pytest.approx(0.0, abs=1e-15)
    assert r["land melt"] == 0.0 and r["land ligroundf"] == 0.0
    assert r["left out"] == pytest.approx(0.0, abs=1e-15)
    assert cmb.failures([{"year": 2300, **r}], tol=1e-6) == []


def test_a_forward_split_that_does_not_reproduce_fails():
    r"""Half of cell 1's front melt written back into libmassbffl: the melt
    is the same, and the inflow's share of it is not what the file says."""
    prev, cur = _year("now")
    half = cur["lifmassbf"][1] / 2
    cur["lifmassbf"][1] -= half
    cur["libmassbffl"][1] += half
    r = cmb.year_booking(prev, cur, W)
    assert r["split"] == pytest.approx(0.5)
    assert [k for _, k, _ in cmb.failures([{"year": 2300, **r}], tol=1e-6)] == ["split"]


def test_a_missing_series_is_an_input_error(tmp_path):
    assert cmb.main([str(tmp_path / "nothing_ismip7_annual.h5")]) == 2
