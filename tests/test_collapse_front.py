r"""``ISMIP7_FRACTURE=mask_front``: a shelf collapses from its front.

The ISMIP7 collapse masks flag the Ross and Filchner-Ronne shelves near their
grounding lines decades before they reach the fronts. Emptying every flagged
floating cell (``mask``) opens holes hundreds of kilometres behind a standing
front, the momentum balance treats a hole as open ocean, and the groups doing
it report 40 % to 100 % more sea level by 2300 from that alone (discussion
#30, September 2026). ``mask_front`` is the other end-member they settled on:
a flagged cell goes once open water has reached it through flagged cells.

The rule is a pure function over the cell arrays, exercised here on a strip
of cells, and then through the real facet-neighbour operator on a small
Firedrake mesh, where "shares a facet" has to mean an edge and not a corner.

The group has to pick one mode from runs under each (issue #10), so every run
says which mode it ran and how many flagged floating cells the mode removed
and held: in the log, in the timeseries and in the committed run record.
"""
import os
import sys

import numpy as np
import pytest

from icepack2_tools.front import (
    COLLAPSE_CSV_COLUMNS, COLLAPSE_MARKER, collapse_banner,
    collapse_cell_counts, collapse_csv_fields, front_connected,
)
from icepack2_tools.runconfig import FRACTURE_MASK_MODES, FRACTURE_MODES, fracture

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _strip_neighbours(mask):
    r"""Cells of a 1-D strip with a neighbour in ``mask``."""
    out = np.zeros_like(mask)
    out[1:] |= mask[:-1]
    out[:-1] |= mask[1:]
    return out


def _cells(text):
    r"""``~`` open water, ``#`` flagged shelf, ``.`` unflagged shelf."""
    return np.array([c == "#" for c in text]), np.array([c == "~" for c in text])


def _reached(text):
    flagged, water = _cells(text)
    got = front_connected(flagged, water, _strip_neighbours)
    return "".join("x" if r else c for c, r in zip(text, got))


def test_a_patch_behind_a_standing_front_is_held():
    assert _reached("~~...###..") == "~~...###.."


def test_a_patch_the_water_touches_goes_whole():
    assert _reached("~~####....") == "~~xxxx...."


def test_an_unflagged_cell_stops_the_collapse():
    # the front-side patch goes; the one past the unflagged cell waits
    assert _reached("~##.###") == "~xx.###"


def test_water_on_both_sides_reaches_from_both():
    assert _reached("~##..##~") == "~xx..xx~"


def test_a_cell_emptied_earlier_is_reached_by_definition():
    r"""After an advance an emptied cell is flagged AND open water. It has to
    stay in the set: that is what names it to the thickness floor, and a
    floored cell would be re-calved on every advance for the rest of the run."""
    flagged = np.array([False, True, True, False])
    water = np.array([True, True, False, False])
    assert front_connected(flagged, water, _strip_neighbours).tolist() == [False, True, True, False]


def test_the_inputs_are_left_alone():
    flagged, water = _cells("~###")
    before = flagged.copy(), water.copy()
    front_connected(flagged, water, _strip_neighbours)
    assert (flagged == before[0]).all() and (water == before[1]).all()


def test_the_sweep_ends_on_the_collective_answer():
    r"""``any_rank`` is the collective OR. A rank with nothing new must keep
    sweeping while another still grows, or the ranks stop calling the
    collective assembly in step and the run hangs."""
    calls = []

    def any_rank(local):
        calls.append(bool(local))
        return len(calls) < 3                     # "another rank" grew twice more

    flagged, water = _cells("...")
    front_connected(flagged, water, _strip_neighbours, any_rank)
    assert calls == [False, False, False]


def test_both_mask_modes_are_mask_modes(monkeypatch):
    assert FRACTURE_MASK_MODES == ("mask", "mask_front")
    monkeypatch.setenv("ISMIP7_FRACTURE", "mask_front")
    assert fracture() == "mask_front"
    monkeypatch.setenv("ISMIP7_FRACTURE", "front")
    with pytest.raises(ValueError, match="mask_front"):
        fracture()


def test_the_control_and_ocx_refuse_either_mask_mode(monkeypatch):
    from icepack2_tools.forcing import reject_collapse_mask
    for mode in FRACTURE_MASK_MODES:
        monkeypatch.setenv("ISMIP7_FRACTURE", mode)
        with pytest.raises(ValueError, match=f"ISMIP7_FRACTURE={mode}"):
            reject_collapse_mask("the control")
    monkeypatch.setenv("ISMIP7_FRACTURE", "none")
    reject_collapse_mask("the control")


# --- what a run reports -------------------------------------------------------

def test_every_mode_prints_its_banner():
    r"""``none`` included: a log with no banner is then a run that predates
    it, and the mode of any later run can be read off its log."""
    for mode in FRACTURE_MODES:
        assert collapse_banner(mode).startswith(f"{COLLAPSE_MARKER} ISMIP7_FRACTURE={mode} (")
    assert "no cell is removed" in collapse_banner("none")


def test_mask_holds_nothing_and_mask_front_holds_the_hole():
    flagged, water = _cells("~##..###..")
    assert collapse_cell_counts(flagged, flagged) == (5, 5, 0)                # mask
    reached = front_connected(flagged, water, _strip_neighbours)
    assert collapse_cell_counts(flagged, reached) == (5, 2, 3)                # mask_front


def test_a_cell_emptied_earlier_still_counts_as_removed():
    r"""It stays flagged and afloat, so it stays in both sets: ``removed``
    accumulates, and flagged = removed + held on every advance."""
    flagged = np.array([False, True, True, True])
    water = np.array([True, True, False, False])                  # cell 1 went last advance
    reached = front_connected(flagged, water, _strip_neighbours)
    assert collapse_cell_counts(flagged, reached) == (3, 3, 0)


def test_every_count_goes_through_the_reduction():
    r"""``.sum()`` on a cell array is one rank's slice (AGENTS.md section 3).
    Each of the three counts is reduced, the held one too."""
    flagged, water = _cells("~##..###..")
    reached = front_connected(flagged, water, _strip_neighbours)

    def two_ranks(mask):                                          # a twin rank holding the same cells
        return 2 * np.count_nonzero(mask)

    assert collapse_cell_counts(flagged, reached, two_ranks) == (10, 4, 6)


def test_without_a_mask_nothing_is_counted_or_reduced():
    r"""Under ``none`` no rank enters the collective, and the counts are zero."""
    def never(mask):
        raise AssertionError("a reduction ran with no mask")

    assert collapse_cell_counts(None, None, never) == (0, 0, 0)


def test_the_timeseries_carries_the_counts_and_an_older_series_keeps_its_width():
    new = "year,vaf_mm_sle,resid_gt,amb_gtyr," + ",".join(COLLAPSE_CSV_COLUMNS) + "\n"
    assert COLLAPSE_CSV_COLUMNS[-1] == "collapse_held_cells"
    assert collapse_csv_fields(new, (5, 2, 3)) == ",5,2,3"
    # a series begun before the columns existed is resumed under its own
    # header, so its rows must not grow past it
    assert collapse_csv_fields("year,vaf_mm_sle,resid_gt,amb_gtyr\n", (5, 2, 3)) == ""


def test_the_run_record_states_the_mode_and_lifts_the_counts(tmp_path, monkeypatch):
    sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))
    import core_report

    monkeypatch.delenv("ISMIP7_FRACTURE", raising=False)
    assert core_report.effective_env()["ISMIP7_FRACTURE"].startswith("none    # default")
    monkeypatch.setenv("ISMIP7_FRACTURE", "mask_front")
    assert core_report.effective_env()["ISMIP7_FRACTURE"] == "mask_front"

    banner = collapse_banner("mask_front")
    closing = f"{COLLAPSE_MARKER} ISMIP7_FRACTURE=mask_front ended t=2301.0 with flagged=5 removed=2 held=3 cells"
    log = tmp_path / "run.log"
    log.write_text(f"  {banner}\n  t=2015.1\n      collapse [cells]: flagged=5 removed=2 held=3\n{closing}\n")
    assert core_report.collapse_record(str(log)) == [banner, closing]
    quiet = tmp_path / "old.log"
    quiet.write_text("  t=2015.1\n")
    (none,) = core_report.collapse_record(str(quiet))
    assert "predates the collapse banner" in none


# --- through the real facet-neighbour operator ------------------------------

@pytest.fixture(scope="module")
def shelf():
    fd = pytest.importorskip("firedrake")
    from icepack2_tools.front import facet_neighbours
    mesh = fd.RectangleMesh(8, 4, 8.0, 4.0, quadrilateral=True)     # 1 x 1 cells
    Q = fd.FunctionSpace(mesh, "DG", 0)
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro
    return facet_neighbours(Q), xy[:, 0], xy[:, 1]


def test_a_neighbour_shares_an_edge_not_a_corner(shelf):
    neighbours_of, x, y = shelf
    one = (np.abs(x - 3.5) < 0.1) & (np.abs(y - 1.5) < 0.1)
    assert one.sum() == 1
    got = neighbours_of(one)
    assert got.sum() == 4 and not got[one].any()
    assert np.allclose(np.abs(x[got] - 3.5) + np.abs(y[got] - 1.5), 1.0)


def test_a_hole_is_held_and_a_corridor_to_the_water_lets_it_go(shelf):
    neighbours_of, x, y = shelf
    water = x < 1.0                                                # the column at the front
    patch = (x > 4.0) & (x < 7.0) & (y > 1.0) & (y < 3.0)          # 3 x 2 cells, well inside
    assert not front_connected(patch, water, neighbours_of).any()
    corridor = (x > 1.0) & (x < 4.0) & (np.abs(y - 1.5) < 0.1)     # flagged cells out to the front
    got = front_connected(patch | corridor, water, neighbours_of)
    assert (got == (patch | corridor)).all()
    # a corridor that only touches the patch at a corner does not connect it
    diagonal = (x > 1.0) & (x < 4.0) & (np.abs(y - 0.5) < 0.1)
    got = front_connected(patch | diagonal, water, neighbours_of)
    assert (got == diagonal).all()


def test_the_held_hole_is_counted_over_the_mesh_communicator(shelf):
    from firedrake import COMM_WORLD
    from icepack2_tools.mpi_stats import global_count
    neighbours_of, x, y = shelf
    water = x < 1.0
    patch = (x > 4.0) & (x < 7.0) & (y > 1.0) & (y < 3.0)          # the 3 x 2 hole
    corridor = (x > 1.0) & (x < 4.0) & (np.abs(y - 1.5) < 0.1)

    def count(mask):
        return global_count(mask, COMM_WORLD)

    held = front_connected(patch, water, neighbours_of)
    assert collapse_cell_counts(patch, held, count) == (6, 0, 6)
    gone = front_connected(patch | corridor, water, neighbours_of)
    assert collapse_cell_counts(patch | corridor, gone, count) == (9, 9, 0)
