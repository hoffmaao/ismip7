r"""The ISMIP7 writer's regridding policies and time encoding, on tiny
synthetic operators (no mesh, no data).

The time stamps are the ones the discussion board settled on (#16, #20):
a state written for year 2015 is stamped 2016-01-01, day 60630 since
1850-01-01 on the standard calendar; a flux for 2015 is stamped 2015-07-01,
day 60446. The fill policies follow the request's csv: ``forbidden``
means the uncovered part of a pixel counts as zero (sums are conserved),
``outside_domain`` means the mean over the covered part, ``no_ice`` and
friends mean over the masked part.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "antarctica", "scripts"))
scipy_sparse = pytest.importorskip("scipy.sparse")
wio = pytest.importorskip("write_ismip7_output")


def _operator():
    # two pixels, three cells: pixel 0 fully covered by cells 0 and 1 (half
    # each), pixel 1 one quarter covered by cell 2
    A = wio.PIXEL_AREA
    return scipy_sparse.csr_matrix(np.array([[A / 2, A / 2, 0.0], [0.0, 0.0, A / 4]]))


def test_forbidden_policy_conserves_the_sum():
    W = _operator(); v = np.array([100.0, 300.0, 400.0])
    out = wio.regrid(W, v, "forbidden")
    assert np.allclose(out, [200.0, 100.0])                   # pixel 1: a quarter of 400
    assert np.isclose((out * wio.PIXEL_AREA).sum(), (W @ v).sum())


def test_outside_domain_policy_is_the_covered_mean():
    W = _operator(); v = np.array([100.0, 300.0, 400.0])
    out = wio.regrid(W, v, "outside_domain")
    assert np.allclose(out, [200.0, 400.0])                   # partial coverage does not dilute an elevation


def test_masked_policy_averages_over_the_mask():
    W = _operator(); v = np.array([10.0, 30.0, 50.0]); mask = np.array([True, False, True])
    out = wio.regrid(W, v, "no_ice", mask)
    assert np.allclose(out, [10.0, 50.0])
    assert np.isnan(wio.regrid(W, v, "no_ice", np.zeros(3, bool))).all()


def test_time_encoding_matches_the_board():
    assert wio.days_since_1850(2016, 1, 1) == 60630           # state written for year 2015
    assert wio.days_since_1850(2015, 7, 1) == 60446           # flux for year 2015
