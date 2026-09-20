r"""The timing matrix leads with what a projection costs, not what a lane took.

A lane times a fixed number of steps at its mesh's own timestep, so lanes are
only comparable per simulated year: seconds_per_step / dt. The 285-year column
is that figure extrapolated over 2015-2300, in hours up to two days and in days
beyond.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "antarctica" / "scripts"))

from build_timing_matrix import (  # noqa: E402
    _full_simulation,
    _minutes_per_year,
    _solver_work,
    _timing_table,
)


def test_a_year_costs_the_step_time_over_the_lane_s_own_dt():
    assert _minutes_per_year({"seconds_per_step": 6.0, "dt": 0.125}) == "0.8"
    assert _minutes_per_year({"seconds_per_step": 6.0, "dt": 0.05}) == "2.0"


@pytest.mark.parametrize("seconds_per_step, expected", [
    (7.5, "4.8 h"),      # 60 s/yr
    (75.0, "47.5 h"),
    (76.0, "2.0 d"),     # past two days the unit changes
    (750.0, "19.8 d"),
])
def test_the_full_simulation_is_285_years_in_hours_or_days(seconds_per_step, expected):
    assert _full_simulation({"seconds_per_step": seconds_per_step, "dt": 0.125}) == expected


def test_a_record_without_a_dt_is_left_blank_rather_than_guessed():
    legacy = {"seconds_per_step": 6.0}
    assert _minutes_per_year(legacy) == "—"
    assert _full_simulation(legacy) == "—"


def test_the_unit_is_in_the_header_only_when_every_cell_shares_it():
    record = {"seconds_per_step": 6.0, "dt": 0.125, "vertices": 10, "cells": 20}
    index = {(2500, 25000, 16): record}
    shared = _timing_table([(2500, 25000)], index, _minutes_per_year, "min/yr")
    mixed = _timing_table([(2500, 25000)], index, _full_simulation)
    assert "16 cores (min/yr)" in shared[0]
    assert "| 16 cores |" in mixed[0]
    assert shared[2].startswith("| 2500 | 25000 | 10 | 20 | 0.8 | — |")


def test_solver_work_is_newton_per_step_by_krylov_per_newton():
    r"""Under GAMG a step's cost is its Krylov count; the matrix shows it beside
    the direct solve's one-per-Newton so the two campaigns can be read together."""
    def summary(newton, krylov):
        return {"diagnostic_solve_summary": {
            "count": 10, "snes_iterations_total": newton, "linear_iterations_total": krylov}}
    # Records from before SCPC counted its own solves hold SNES's count only.
    assert _solver_work(summary(42, 42)) == "4.2 × 1.0†"
    assert _solver_work(summary(42, 1512)) == "4.2 × 36.0†"
    assert _solver_work(summary(0, 0)) == "—"
    assert _solver_work({"diagnostic_solve_summary": {"reason_counts": {"2": 10}}}) == "—"
    assert _solver_work({}) == "—"
