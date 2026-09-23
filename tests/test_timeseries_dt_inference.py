r"""The audit divides per-step budget columns by the timestep, so inferring
the timestep wrongly rescales every rate it reports. A dt=0.05 run whose year
column is written to one decimal repeats years, and the median difference then
lands on exactly twice the truth."""
import numpy as np
import pytest

from antarctica.scripts.check_ismip6_track import infer_dt


def test_a_uniform_series_gives_its_own_step():
    yr = 2015.0 + 0.1 * np.arange(40)
    assert infer_dt(yr) == pytest.approx(0.1)


def test_a_rounded_column_does_not_double_the_step():
    r"""The defect: 40 steps of 0.05 yr written to one decimal produce
    differences of 0.0 and 0.1, whose median is 0.1."""
    exact = 2015.0 + 0.05 * np.arange(40)
    rounded = np.round(exact, 1)
    assert float(np.median(np.diff(rounded))) == pytest.approx(0.1)
    # the fallback is the mean spacing, which rounding perturbs by O(1/n);
    # what matters is that it is near 0.05 rather than a factor of two out
    assert infer_dt(rounded) == pytest.approx(0.05, rel=0.05)


def test_a_finer_step_survives_the_same_rounding():
    exact = 2015.0 + 0.025 * np.arange(81)
    assert infer_dt(np.round(exact, 1)) == pytest.approx(0.025, rel=0.05)


def test_a_single_row_falls_back_to_a_year():
    assert infer_dt(np.array([2015.0])) == pytest.approx(1.0)


def test_the_year_column_is_written_finely_enough_to_resolve_the_step():
    r"""Fixing the reader is not enough: the writer has to stop discarding the
    information. Guard the format string itself, since nothing else in the
    suite runs a sub-0.1 transient to disk."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1]
           / "antarctica" / "scripts" / "simulation.py").read_text()
    assert 'f"{row[0]:.4f},' in src, "the timeseries year column lost precision"
    assert 'f"{row[0]:.1f},' not in src
