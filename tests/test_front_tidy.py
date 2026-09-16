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

from icepack2_tools.front import clamp_thickness, retreat_slivers

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


# ---- the h_clamp floor ---------------------------------------------------
#
# Exercises :func:`icepack2_tools.front.clamp_thickness`. Every rule that
# empties a cell has to name those cells as ice-free here; miss one and the
# floor refills the cell out of nothing on the next advance and the emptying
# rule removes it again, fabricating h_clamp per step of clamp mass and of
# the calving flux the submission reports.

CLAMP = 10.0


def test_the_floor_lifts_a_thin_ice_cell():
    h = np.array([0.0, 3.0, 40.0])
    clamp_thickness(h, CLAMP)
    assert np.allclose(h, [CLAMP, CLAMP, 40.0])


def test_an_ice_free_cell_keeps_its_zero_floor():
    h = np.array([0.0, 0.0])
    clamp_thickness(h, CLAMP, np.array([True, False]))
    assert np.allclose(h, [0.0, CLAMP])


def test_a_collapsed_shelf_cell_stays_empty_across_advances():
    r"""The ISMIP7 collapse mask empties floating cells. Without the
    exemption the floor refills each one every advance and the collapse
    removal re-calves it, so licalvf carries CLAMP metres per step forever."""
    collapse = np.array([True, True, False])
    grounded = np.array([False, True, False])
    collapsed = collapse & ~grounded          # floating and flagged: cell 0
    recalved = 0.0
    h = np.array([0.0, 0.0, 0.0])
    for _ in range(5):
        clamp_thickness(h, CLAMP, None, None, collapsed)
        hit = collapsed & (h > 0.0)
        recalved += float(h[hit].sum())
        h[hit] = 0.0
    assert recalved == 0.0                    # nothing fabricated to re-calve
    assert h[0] == 0.0                        # the collapsed cell stays empty
    assert h[1] == CLAMP                      # grounded: the floor still applies
    assert h[2] == CLAMP                      # untouched by the mask


def test_the_exempt_masks_are_a_union():
    h = np.zeros(4)
    clamp_thickness(h, CLAMP,
                    np.array([True, False, False, False]),
                    np.array([False, True, False, False]),
                    np.array([False, False, True, False]))
    assert np.allclose(h, [0.0, 0.0, 0.0, CLAMP])
