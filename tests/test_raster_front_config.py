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


def test_the_front_sampling_names_its_mask_year_and_samples_rasters_by_vertex():
    assert "greene2015" in R.RASTER_SAMPLES
    assert R.raster_front_year("greene2015") == 2015
    assert R.raster_front_year("vertex") is None
    assert R.raster_front_year("cell_mean") is None
    assert R.raster_base_method("greene2015") == "vertex"
    assert R.raster_base_method("cell_mean") == "cell_mean"
    assert R.target_mesh_geometry_method("vertex") == R.TARGET_MESH_GEOMETRY_METHOD
    assert R.target_mesh_geometry_method("greene2015") != R.TARGET_MESH_GEOMETRY_METHOD
    assert "2015" in R.target_mesh_geometry_method("greene2015")


def test_a_forward_on_its_maps_mesh_follows_the_map(monkeypatch):
    monkeypatch.delenv("ISMIP7_RASTER_SAMPLE", raising=False)
    assert R.forward_raster_sample("greene2015", transfer=False) == "greene2015"
    assert R.forward_raster_sample(b"vertex", transfer=False) == "vertex"
    assert R.forward_raster_sample(None, transfer=False) == "vertex"
    monkeypatch.setenv("ISMIP7_RASTER_SAMPLE", "greene2015")
    assert R.forward_raster_sample("greene2015", transfer=False) == "greene2015"
    with pytest.raises(RuntimeError, match="follows the MAP"):
        R.forward_raster_sample("vertex", transfer=False)


def test_a_transfer_rebuilds_with_this_runs_sampling(monkeypatch):
    monkeypatch.delenv("ISMIP7_RASTER_SAMPLE", raising=False)
    assert R.forward_raster_sample("vertex", transfer=True) == R.RASTER_SAMPLE_DEFAULT
    monkeypatch.setenv("ISMIP7_RASTER_SAMPLE", "greene2015")
    assert R.forward_raster_sample("vertex", transfer=True) == "greene2015"


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
                                {"raster_sample": "greene2015"})
    # a warm start older than the record is no mismatch
    assert not objective_mismatches({}, {"raster_sample": "greene2015"})
