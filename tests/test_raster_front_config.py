"""Which raster sampling a run builds its geometry with, and which melt
calibration it takes for it (issue #167).

A forward on its MAP's mesh uses the MAP's geometry, so the sampling the MAP
records stands. A forward on another mesh rebuilds the geometry from
BedMachine there, with this run's sampling. A run melts with the tracked
calibration fitted under its own sampling, and a sampling without one is
refused.
"""

import pytest

from icepack2_tools import runconfig as R
from icepack2_tools.handoff import OBJECTIVE_KEYS, objective_mismatches


def test_the_front_sampling_samples_rasters_by_vertex():
    assert "vertex_front" in R.RASTER_SAMPLES
    assert R.raster_front("vertex_front")
    assert not R.raster_front("vertex") and not R.raster_front("cell_mean")
    assert R.raster_base_method("vertex_front") == "vertex"
    assert R.raster_base_method("cell_mean") == "cell_mean"
    assert R.target_mesh_geometry_method("vertex") == R.TARGET_MESH_GEOMETRY_METHOD
    assert R.target_mesh_geometry_method("vertex_front") != R.TARGET_MESH_GEOMETRY_METHOD


def test_a_forward_on_its_maps_mesh_follows_the_map(monkeypatch):
    monkeypatch.delenv("ISMIP7_RASTER_SAMPLE", raising=False)
    assert R.forward_raster_sample("vertex_front", transfer=False) == "vertex_front"
    assert R.forward_raster_sample(b"vertex", transfer=False) == "vertex"
    assert R.forward_raster_sample(None, transfer=False) == "vertex"
    monkeypatch.setenv("ISMIP7_RASTER_SAMPLE", "vertex_front")
    assert R.forward_raster_sample("vertex_front", transfer=False) == "vertex_front"
    with pytest.raises(RuntimeError, match="follows the MAP"):
        R.forward_raster_sample("vertex", transfer=False)


def test_a_transfer_rebuilds_with_this_runs_sampling(monkeypatch):
    monkeypatch.delenv("ISMIP7_RASTER_SAMPLE", raising=False)
    assert R.forward_raster_sample("vertex", transfer=True) == R.RASTER_SAMPLE_DEFAULT
    monkeypatch.setenv("ISMIP7_RASTER_SAMPLE", "vertex_front")
    assert R.forward_raster_sample("vertex", transfer=True) == "vertex_front"


def test_a_run_takes_the_calibration_of_its_sampling(monkeypatch):
    for knob in ("ISMIP7_DELTAT_PER_BASIN_NPZ", "ISMIP7_K_PER_BASIN_NPZ",
                 "ISMIP7_RASTER_SAMPLE", "ISMIP7_K_SCALE"):
        monkeypatch.delenv(knob, raising=False)
    assert R.deltat_per_basin_npz("vertex") == R.MELT_CALIBRATIONS["vertex"]
    assert R.deltat_per_basin_npz() == R.MELT_CALIBRATION_DEFAULT
    missing = [m for m in R.RASTER_SAMPLES if m not in R.MELT_CALIBRATIONS]
    for method in missing:
        with pytest.raises(FileNotFoundError, match=method):
            R.deltat_per_basin_npz(method)


def test_the_sampling_is_part_of_the_objective():
    assert "raster_sample" in OBJECTIVE_KEYS
    assert objective_mismatches({"raster_sample": "vertex"},
                                {"raster_sample": "vertex_front"})
    # a warm start older than the record is no mismatch
    assert not objective_mismatches({}, {"raster_sample": "vertex_front"})


def test_a_warm_start_state_comes_with_its_geometry_or_as_a_first_guess(monkeypatch):
    r"""Issue #167: a warm start sampled another way keeps its geometry, and
    with it by default its mixed state; ISMIP7_WARM_START_STATE=1 loads the
    state as the first guess on the warm start's own mesh only."""
    monkeypatch.delenv("ISMIP7_WARM_START_STATE", raising=False)
    assert R.warm_start_state(geometry_taken=True, same_mesh=True) == (True, False)
    assert R.warm_start_state(geometry_taken=False, same_mesh=True) == (False, False)
    monkeypatch.setenv("ISMIP7_WARM_START_STATE", "1")
    assert R.warm_start_state(geometry_taken=False, same_mesh=True) == (True, True)
    assert R.warm_start_state(geometry_taken=True, same_mesh=True) == (True, False)
    with pytest.raises(ValueError, match="own mesh"):
        R.warm_start_state(geometry_taken=False, same_mesh=False)
    monkeypatch.setenv("ISMIP7_WARM_START_STATE", "0")
    assert R.warm_start_state(geometry_taken=False, same_mesh=False) == (False, False)
