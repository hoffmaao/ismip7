r"""check_ismip6_track's runaway detector: sustained growth or a sustained
high year, never a single-step spike (icepack/ismip7#33)."""
import importlib.util
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _detect():
    spec = importlib.util.spec_from_file_location(
        "chk", os.path.join(REPO, "antarctica", "scripts", "check_ismip6_track.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.runaway_detected


def test_a_steady_run_passes():
    assert _detect()(np.full(300, 1000.0), 0.1) is False


def test_one_emptying_step_is_not_a_runaway():
    """An emptying event books one step of 60,000 Gt/yr and the budget closes
    (core 2 of the July matrix); the old peak clause failed it."""
    d = np.full(300, 1000.0)
    d[150] = 60000.0
    assert _detect()(d, 0.1) is False


def test_one_fast_year_that_settles_is_not_a_runaway():
    """Growth of 2.4x into one year and back (core 3): sustained means two."""
    d = np.full(300, 1000.0)
    d[100:110] = 2400.0
    assert _detect()(d, 0.1) is False


def test_the_july_blow_up_still_fails_in_year_one():
    """~2x per 0.1-yr step from 1000 Gt/yr: the median of the first year is
    already above 6000 and the growth is sustained."""
    d = 1000.0 * 2.0 ** np.arange(30)
    assert _detect()(d, 0.1) is True


def test_two_consecutive_years_of_growth_fail_even_below_the_ceiling():
    d = np.concatenate([np.full(10, 1000.0), np.full(10, 1600.0), np.full(10, 2600.0), np.full(10, 2600.0)])
    assert _detect()(d, 0.1) is True


def test_a_ten_step_probe_that_sits_above_the_ceiling_fails():
    """The 32 km probes that blew up on step 1 have one year block; its median
    carries the verdict."""
    d = np.array([1200.0, 30000.0, 40000.0, 45000.0, 48000.0, 50000.0, 52000.0, 55000.0, 58000.0, 60000.0])
    assert _detect()(d, 0.1) is True


def test_dt_does_not_change_the_answer():
    """Year blocks: the same spike at dt=0.05 is still one step of one year."""
    d = np.full(600, 1000.0)
    d[300] = 60000.0
    assert _detect()(d, 0.05) is False


def test_a_spike_on_a_one_step_trailing_block_is_not_a_runaway():
    """A t=0 row leaves a one-step final block; its median would be the spike."""
    d = np.full(301, 1000.0)
    d[-1] = 60000.0
    assert _detect()(d, 0.1) is False


def test_a_full_trailing_year_above_the_ceiling_still_fails():
    d = np.concatenate([np.full(290, 1000.0), np.full(10, 8000.0)])
    assert _detect()(d, 0.1) is True
