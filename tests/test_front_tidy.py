r"""The front tidy rule: what a transport advance leaves in a cell.

Exercises :func:`icepack2_tools.front.retreat_slivers`, the rule that separates
a retreating front's remainder (removed, tallied as calving) from the
sub-threshold inflow that carries a free front's ADVANCE (never touched).
Drop the ``h_old > front_hmin`` conjunct and a configured calving law is
silently pinned at its starting extent, which is what these tests catch.

Serial, no data files, no run setup.
"""

import numpy as np
import pytest

from icepack2_tools.front import retreat_slivers

HMIN = 1.0


def _one(h_new, h_old, front_hmin=HMIN):
    r"""The rule applied to a single cell, as a bool."""
    return bool(retreat_slivers(np.array([h_new]), np.array([h_old]),
                                front_hmin)[0])


def test_subthreshold_inflow_survives():
    r"""An ice-free cell that just received 0.05 m keeps it: this is advance."""
    assert not _one(h_new=0.05, h_old=0.0)


def test_inflow_accumulates_until_it_joins_the_extent():
    r"""Repeated sub-threshold inflow crosses the threshold instead of being
    destroyed each step, so the front can advance into the cell."""
    h = 0.0
    for _ in range(20):
        h_old, h = h, h + 0.05
        if _one(h_new=h, h_old=h_old):
            h = 0.0
    assert h == pytest.approx(1.0)
    assert not _one(h_new=h, h_old=0.95)


def test_retreat_remainder_removed():
    r"""A cell that held ice and was shed down to a sliver is taken."""
    assert _one(h_new=0.5, h_old=5.0)


def test_ice_above_threshold_untouched():
    assert not _one(h_new=3.0, h_old=5.0)


def test_emptied_cell_untouched():
    r"""Already at zero: nothing to remove, so it must not be selected."""
    assert not _one(h_new=0.0, h_old=5.0)


@pytest.mark.parametrize("h_new,h_old,expected", [
    (1.0, 5.0, True),     # exactly at the threshold, was ice -> a sliver
    (1.0, 0.0, False),    # exactly at the threshold, was ice-free -> advance
    (0.5, 1.0, False),    # h_old exactly at the threshold is not "was ice"
    (1.0001, 5.0, False),  # just above the threshold -> still ice
])
def test_threshold_boundaries(h_new, h_old, expected):
    assert _one(h_new=h_new, h_old=h_old) is expected


def test_mask_over_many_cells():
    r"""Vectorised over a cell array, the shape the transport passes in."""
    h_old = np.array([0.0, 0.0, 5.0, 5.0, 5.0, 0.5])
    h_new = np.array([0.05, 0.0, 0.5, 3.0, 0.0, 0.7])
    got = retreat_slivers(h_new, h_old, HMIN)
    np.testing.assert_array_equal(
        got, np.array([False, False, True, False, False, False]))
    assert got.dtype == np.bool_


def test_removal_conserves_the_tally():
    r"""What the mask removes is exactly what gets booked to calving."""
    area = np.array([2.0, 2.0, 3.0, 3.0])
    h_old = np.array([0.0, 5.0, 5.0, 5.0])
    h = np.array([0.05, 0.5, 4.0, 0.25])
    before = float((h * area).sum())
    sliver = retreat_slivers(h, h_old, HMIN)
    calved = float((h[sliver] * area[sliver]).sum())
    h[sliver] = 0.0
    assert calved == pytest.approx(0.5 * 2.0 + 0.25 * 3.0)
    assert float((h * area).sum()) == pytest.approx(before - calved)


def test_on_a_dg0_unit_mesh():
    r"""The same rule over real DG0 cell data, the arrays `_advance` passes."""
    fd = pytest.importorskip("firedrake")

    mesh = fd.UnitSquareMesh(2, 2)
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    h_new, h_old = fd.Function(Q0), fd.Function(Q0)
    n = h_new.dat.data.shape[0]
    assert n >= 4

    start, post = np.zeros(n), np.zeros(n)
    start[0], post[0] = 100.0, 90.0   # ice, still thick
    start[1], post[1] = 5.0, 0.5      # ice, shed to a sliver
    start[2], post[2] = 0.0, 0.05     # ice-free, inflow
    h_old.dat.data[:] = start
    h_new.dat.data[:] = post

    area = fd.assemble(fd.TestFunction(Q0) * fd.dx).dat.data_ro.copy()
    m_before = float(fd.assemble(h_new * fd.dx))

    data = h_new.dat.data
    sliver = retreat_slivers(data, h_old.dat.data_ro, HMIN)
    calved = float((data[sliver] * area[sliver]).sum())
    data[sliver] = 0.0

    after = h_new.dat.data_ro
    assert after[2] == pytest.approx(0.05)
    assert after[1] == 0.0
    assert after[0] == pytest.approx(90.0)
    assert float(fd.assemble(h_new * fd.dx)) == pytest.approx(m_before - calved)
