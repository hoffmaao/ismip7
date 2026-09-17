r"""The roots a second checkout does not carry.

A cluster can hold more than one checkout of this repository. The code moves
with the invocation; the forcing tree, the meshes, the MAPs, the observational
rasters and the calibration npz files do not - they are gitignored, so a fresh
clone has none of them. These check the Python half of that split (site_env.sh
and tests/test_site_core.py cover the shell half), and in particular that every
default still resolves into the invoking checkout when ISMIP7_SHARE is unset,
so a site holding a single checkout is provably unaffected.
"""
import os

import pytest

from icepack2_tools.runconfig import (
    k_per_basin_candidates, obs_data_root, shared_results_root,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANT = os.path.join(REPO, "antarctica")


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for key in ("ISMIP7_SHARE", "ISMIP7_OBS_DATA_ROOT", "ISMIP7_K_PER_BASIN_NPZ"):
        monkeypatch.delenv(key, raising=False)


def test_the_obs_root_is_this_checkout_when_nothing_says_otherwise():
    assert obs_data_root() == os.path.join(ANT, "data")


def test_the_obs_root_follows_its_variable(monkeypatch):
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", "/projects/shared/antarctica/data")
    assert obs_data_root() == "/projects/shared/antarctica/data"


def test_there_is_no_shared_results_root_without_a_share():
    assert shared_results_root() is None


def test_the_shared_results_root_hangs_off_the_share(monkeypatch):
    monkeypatch.setenv("ISMIP7_SHARE", "/projects/ah301/ismip7")
    assert shared_results_root() == "/projects/ah301/ismip7/antarctica/results"


def test_k_is_looked_for_in_this_checkout_first_then_the_share(monkeypatch):
    """Order matters: a checkout that has its own calibration must win over the
    shared one, or a deliberate local recalibration is silently ignored."""
    monkeypatch.setenv("ISMIP7_SHARE", "/shared")
    got = k_per_basin_candidates("/here/results", 2000)
    assert got == [
        "/here/results/calibrated_K_per_basin_2000.npz",
        "/here/results/calibrated_K_per_basin_2500.npz",
        "/shared/antarctica/results/calibrated_K_per_basin_2000.npz",
        "/shared/antarctica/results/calibrated_K_per_basin_2500.npz",
    ]


def test_without_a_share_only_this_checkout_is_searched():
    assert k_per_basin_candidates("/here/results", 2000) == [
        "/here/results/calibrated_K_per_basin_2000.npz",
        "/here/results/calibrated_K_per_basin_2500.npz",
    ]


def test_the_share_is_not_searched_twice_when_it_is_this_checkout(monkeypatch):
    """A site whose ISMIP7_SHARE is the checkout repeats the root, and at
    lc=2500 the mesh name and the fallback name coincide. Neither may produce a
    duplicate: the list is what a caller reports when nothing is found."""
    monkeypatch.setenv("ISMIP7_SHARE", REPO)
    results = os.path.join(ANT, "results")
    two_names = k_per_basin_candidates(results, 2000)
    assert two_names == [
        os.path.join(results, "calibrated_K_per_basin_2000.npz"),
        os.path.join(results, "calibrated_K_per_basin_2500.npz"),
    ]
    assert k_per_basin_candidates(results, 2500) == [
        os.path.join(results, "calibrated_K_per_basin_2500.npz")
    ]


def test_an_explicit_npz_is_the_only_candidate(monkeypatch):
    """The override names one file; falling back past it would load a
    calibration the operator did not ask for."""
    monkeypatch.setenv("ISMIP7_SHARE", "/shared")
    monkeypatch.setenv("ISMIP7_K_PER_BASIN_NPZ", "/tmp/mine.npz")
    assert k_per_basin_candidates("/here/results", 2000) == ["/tmp/mine.npz"]


def test_the_inversion_driver_reads_the_obs_root(monkeypatch):
    """It is the driver that did NOT read it: submissions from a second Rice
    checkout died at 'Loading data...' looking for BedMachine while the shell
    layer had named the right tree."""
    import importlib
    import icepack2_tools.runconfig as rcmod
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", "/projects/shared/antarctica/data")
    importlib.reload(rcmod)
    assert rcmod.obs_data_root() == "/projects/shared/antarctica/data"
