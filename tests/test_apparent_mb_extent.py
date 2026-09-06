r"""The frozen apparent-MB reference follows the LIVE ice extent.

Exercises :func:`icepack2_tools.front.clear_reference_where_ice_free`, the rule
that keeps ``a_ref`` off cells that hold no ice. Without it a free calving law
cannot retreat: ``a_ref`` at a t=0 front cell is the terminus outflow, large
and positive, so a cell the front calves is regrown by the transport source on
the very next step and re-booked to calving, oscillating at the t=0 front.

The reference is the DG0 flux divergence, so these are the per-cell arrays the
transport actually adds to its source term.

Serial, no data files, no run setup.
"""

import numpy as np

from icepack2_tools.front import clear_reference_where_ice_free


def test_calved_cell_inside_t0_extent_is_not_regrown():
    r"""The reported failure: a cell INSIDE the t=0 extent that the front has
    just calved keeps a large positive reference, so the next step's source
    regrows it above the extent threshold."""
    front_hmin, dt = 1.0, 0.1
    # Cell 1 sat at the t=0 front, so its reference is the terminus outflow.
    a_ref = np.array([0.5, 40.0, 0.3])
    h = np.array([500.0, 0.0, 480.0])          # cell 1 was just calved
    ice_free = np.array([False, True, False])

    # Before the fix a_ref was masked only at t=0, so cell 1 kept 40 m/yr.
    regrown = h + a_ref * dt
    assert regrown[1] > front_hmin, "guard: the un-masked source does regrow it"

    clear_reference_where_ice_free(a_ref, ice_free)
    assert a_ref[1] == 0.0
    assert (h + a_ref * dt)[1] == 0.0


def test_cells_that_still_hold_ice_keep_their_reference():
    r"""The rule is targeted: it must not disturb the balancing reference on
    the ice, which is what holds a control run steady."""
    a_ref = np.array([0.5, 40.0, -2.0])
    clear_reference_where_ice_free(
        a_ref, np.array([False, True, False]))
    assert a_ref.tolist() == [0.5, 0.0, -2.0]


def test_reference_stays_zero_when_the_front_readvances():
    r"""A cell that re-enters the ice does NOT get its reference back: the
    frozen reference was only ever defined on the t=0 ice."""
    a_ref = np.array([40.0])
    clear_reference_where_ice_free(a_ref, np.array([True]))
    # The front advances back over it: phi turns negative, so it is no longer
    # reported ice-free, and the rule simply does not touch it again.
    clear_reference_where_ice_free(a_ref, np.array([False]))
    assert a_ref[0] == 0.0


def test_rewind_restores_a_reference_the_abandoned_attempt_cleared():
    r"""A rescue-ladder retry must put back the reference on a cell the
    accepted trajectory never calves.

    The step's rewind machinery snapshots the reference at entry and restores
    it alongside h, z and phi. Without that, the dt=0.1 attempt that calved
    cell 1 leaves a_ref[1] = 0 even though the accepted dt=0.025 subcycles
    keep ice there, and the cell then thins at the full unbalanced flux
    divergence for the rest of the run.
    """
    a_ref = np.array([0.5, 40.0, 0.3])
    entry = a_ref.copy()                       # snapshot at step entry

    # The abandoned dt=0.1 attempt: the front passes cell 1, so it is cleared.
    clear_reference_where_ice_free(a_ref, np.array([False, True, False]))
    assert a_ref[1] == 0.0

    # The solve fails, the ladder rewinds and retries at dt/4.
    a_ref[:] = entry
    # At the smaller increment the front does not reach cell 1: nothing clears.
    clear_reference_where_ice_free(a_ref, np.zeros(3, dtype=bool))

    assert a_ref[1] == 40.0, "the retained cell must keep its reference"
    assert a_ref.tolist() == entry.tolist()


def test_nothing_ice_free_is_a_no_op():
    r"""An interior step with a stationary front leaves the field alone."""
    a_ref = np.array([0.5, 40.0, -2.0])
    before = a_ref.copy()
    clear_reference_where_ice_free(a_ref, np.zeros(3, dtype=bool))
    assert np.array_equal(a_ref, before)
