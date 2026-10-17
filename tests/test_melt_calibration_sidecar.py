r"""The offsets fits write the sidecar the forward reads (issue 145).

calibrate_deltaT.py and select_melt_parameters.py write ``<name>.source.json``
beside every offsets file: the npz's sha256, the mesh and its counts, the
raster sampling, every input by name and sha256, the job and the commit. The
forward reads it (``runconfig.melt_calibration_contract``) to refuse a file
that changed, to check the raster sampling and to name the mesh the offsets
were fitted on. Before the fits wrote it, both sidecars were written by hand.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys

import numpy as np
import pytest

try:
    # dual_friction first: it pulls icepack2 -> irksome, which must be
    # imported before any UFL form is assembled.
    import icepack2_tools.dual_friction  # noqa: F401
except ImportError:
    pass
from icepack2_tools.forcing import (                                  # noqa: E402
    check_melt_contract, describe_melt_calibration,
)
from icepack2_tools.melt_selection import TFRule                      # noqa: E402
from icepack2_tools.runconfig import (                                # noqa: E402
    MELT_CALIBRATION_DEFAULT, MELT_CALIBRATION_REQUIRED,
    check_melt_calibration_record, deltat_per_basin_npz, file_sha256,
    melt_calibration_contract, melt_calibration_sidecar,
    refuse_tracked_calibration_out, write_melt_calibration_sidecar,
)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "antarctica", "scripts")
TRACKED = os.path.dirname(MELT_CALIBRATION_DEFAULT)
MESH = "antarctica_250000_25000_buffered20000"


@pytest.fixture
def clean(monkeypatch):
    for k in ("ISMIP7_DELTAT_PER_BASIN_NPZ", "ISMIP7_K_PER_BASIN_NPZ",
              "ISMIP7_K_SCALE", "ISMIP7_K_MELT", "ISMIP7_MELT_SLOPE",
              "ISMIP7_SIN_ALPHA_ANT", "ISMIP7_GEOMETRY_SPACE",
              "ISMIP7_RASTER_SAMPLE", "ISMIP7_SITE", "SLURM_JOB_ID",
              "SLURM_JOB_PARTITION"):
        monkeypatch.delenv(k, raising=False)


@pytest.fixture
def calibrate_melt(monkeypatch, clean):
    r"""``calibrate_melt`` imported afresh, and dropped again afterwards."""
    for name in ("firedrake", "rasterio", "icepack"):
        pytest.importorskip(name)
    monkeypatch.syspath_prepend(SCRIPTS)
    sys.modules.pop("calibrate_melt", None)
    import calibrate_melt
    yield calibrate_melt
    sys.modules.pop("calibrate_melt", None)


def _npz(tmp_path, dT=(0.3, -0.5)):
    path = str(tmp_path / "deltaT_per_basin_25000_K6.500e-05.npz")
    np.savez(path, basin_ids=np.array([3, 9]), deltaT_basin=np.array(dT),
             K=6.5e-5, melt_slope="ant", sin_alpha_ant=5.115e-3,
             geometry_space="dg0")
    return path


def _record(**extra):
    r"""A fit's record, with numpy values where a fit has them."""
    record = {
        "K": 6.5e-5, "site": "iu_quartz", "partition": "debug", "ranks": 1,
        "job": "10669722", "code": "b" * 40, "mesh": MESH,
        "vertices": np.int64(4509), "cells": np.int64(7615),
        "floating_cells": 2923, "geometry_space": "dg0",
        "raster_sample": "vertex", "melt_slope": "ant",
        "sin_alpha_ant": 5.115e-3,
        "obs_table": {"name": "Melt_Paolo_Davison_Adusumilli_imbie2.csv",
                      "sha256": "c" * 64, "total_gtyr": 1067.4},
        "tf_rule": TFRule().as_dict(), "rule_admits": np.bool_(True),
        "melt_total_gtyr": np.float64(1067.386), "unrooted": [np.int64(5)],
    }
    record.update(extra)
    return record


def test_the_sidecar_names_the_npz_as_written(tmp_path):
    npz = _npz(tmp_path)
    path = write_melt_calibration_sidecar(npz, _record())
    assert path == melt_calibration_sidecar(npz)
    with open(path) as f:
        text = f.read()
    assert text.endswith("}\n") and not os.path.exists(path + ".tmp")
    contract = json.loads(text)
    assert list(contract)[:2] == ["file", "sha256"]
    assert contract["file"] == os.path.basename(npz)
    assert contract["sha256"] == file_sha256(npz)
    # numpy values are written as plain JSON, and the rule's pairs as the
    # lists the tracked npz's own tf_rule string reads back as
    assert contract["vertices"] == 4509 and contract["rule_admits"] is True
    assert contract["unrooted"] == [5]
    assert contract["tf_rule"] == json.loads(json.dumps(TFRule().as_dict()))


def test_a_written_sidecar_satisfies_every_reader(clean, monkeypatch, tmp_path):
    npz = _npz(tmp_path)
    write_melt_calibration_sidecar(npz, _record())
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", npz)
    assert deltat_per_basin_npz() == npz
    # the raster sampling is checked, and the mesh named
    assert check_melt_contract(npz, {"raster_sample": "vertex"}) == MESH
    with pytest.raises(ValueError, match="raster_sample=vertex"):
        check_melt_contract(npz, {"raster_sample": "cell_mean"})
    line, = describe_melt_calibration(npz, mesh_basename=MESH + ".msh")
    assert f"fitted on {MESH}; this run's mesh is {MESH}" in line
    assert "does not record" not in line
    # a refit at the tracked K keeps the K's name in the line
    write_melt_calibration_sidecar(npz, _record(
        selected_as="K50", refit_of={"file": "deltaT_per_basin_1000_K6.500e-05.npz",
                                     "sha256": "4" * 64}))
    line, = describe_melt_calibration(npz)
    assert "K 6.500e-05 (K50)" in line
    # the sidecar pins the npz it was written for
    _npz(tmp_path, dT=(0.31, -0.5))
    with pytest.raises(ValueError, match="does not match the sha256"):
        deltat_per_basin_npz()


def test_a_record_the_forward_could_not_trust_is_refused(tmp_path):
    npz = _npz(tmp_path)
    without_sampling = {k: v for k, v in _record().items() if k != "raster_sample"}
    for bad, said in (
            (without_sampling, "lacks raster_sample"),
            (_record(mesh=None), "leaves mesh empty"),
            (_record(written_by="calibrate_deltaT.py --out /N/scratch/someone/refit"),
             "names a path"),
            (_record(mesh_file={"name": "~/meshes/x.msh"}), "names a path"),
            (_record(melt_total_gtyr=float("nan")), "Out of range float"),
            (_record(sin_alpha_cap=float("inf")), "Out of range float")):
        with pytest.raises(ValueError, match=said):
            write_melt_calibration_sidecar(npz, bad)
        assert not os.path.exists(melt_calibration_sidecar(npz))
    # a fit started outside Slurm has no job, and prose may hold a slash
    write_melt_calibration_sidecar(npz, _record(job=None, note="1000 m / 10 km"))
    assert melt_calibration_contract(npz)["job"] is None


def test_the_tracked_sidecar_meets_the_fits_own_schema():
    contract = melt_calibration_contract(MELT_CALIBRATION_DEFAULT)
    assert set(MELT_CALIBRATION_REQUIRED) <= set(contract)
    check_melt_calibration_record(contract)


def test_no_fit_writes_into_the_tracked_calibration(tmp_path):
    before = file_sha256(melt_calibration_sidecar(MELT_CALIBRATION_DEFAULT))
    for directory in (TRACKED, os.path.join(TRACKED, os.pardir, "calibration", "")):
        with pytest.raises(ValueError, match="holds the tracked melt calibration"):
            refuse_tracked_calibration_out(directory)
    refuse_tracked_calibration_out(str(tmp_path))
    with pytest.raises(ValueError, match="holds the tracked melt calibration"):
        write_melt_calibration_sidecar(
            os.path.join(TRACKED, os.path.basename(MELT_CALIBRATION_DEFAULT)),
            _record())
    assert file_sha256(melt_calibration_sidecar(MELT_CALIBRATION_DEFAULT)) == before


def test_both_fits_refuse_the_tracked_directory_before_reading(
        calibrate_melt, monkeypatch):
    def read(*args, **kwargs):
        raise AssertionError("the fit read its inputs")

    spec = importlib.util.spec_from_file_location(
        "calibrate_deltaT", os.path.join(SCRIPTS, "calibrate_deltaT.py"))
    cd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cd)
    monkeypatch.setattr(cd.cm, "_announce_obs_table", read)
    monkeypatch.setattr(sys, "argv", ["calibrate_deltaT.py", "--K", "6.5e-5",
                                      "--out", TRACKED])
    with pytest.raises(ValueError, match="holds the tracked melt calibration"):
        cd.main()
    spec = importlib.util.spec_from_file_location(
        "select_melt_parameters",
        os.path.join(SCRIPTS, "select_melt_parameters.py"))
    sm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sm)
    monkeypatch.setattr(calibrate_melt, "mesh_name", read)
    with pytest.raises(ValueError, match="holds the tracked melt calibration"):
        sm.run_mesh(argparse.Namespace(out=TRACKED))


def test_the_mesh_name_is_the_one_the_forward_uses(calibrate_melt, monkeypatch,
                                                   tmp_path):
    import firedrake as fd
    cm = calibrate_melt
    monkeypatch.setattr(cm, "INV_H5", f"/m/{MESH}.msh")
    assert cm.mesh_name() == MESH
    mesh = fd.UnitSquareMesh(2, 2)

    def checkpoint(name, **attrs):
        path = str(tmp_path / name)
        with fd.CheckpointFile(path, "w") as chk:
            chk.save_mesh(mesh)
            for key, value in attrs.items():
                chk.set_attr("/", key, value)
        return path

    monkeypatch.setattr(cm, "INV_H5", checkpoint("map.h5", mesh_basename=MESH + ".msh"))
    assert cm.mesh_name() == MESH
    # a checkpoint older than the attribute: rebuilt as simulation.py does
    monkeypatch.setattr(cm, "INV_H5", checkpoint(
        "old.h5", lc_coarse=250000, lc=25000, buffer_m=20000.0))
    assert cm.mesh_name() == MESH
    monkeypatch.setattr(cm, "INV_H5", checkpoint("bare.h5"))
    with pytest.raises(ValueError, match="records neither mesh_basename"):
        cm.mesh_name()


def test_the_record_counts_the_mesh_and_names_files_by_basename(
        calibrate_melt, monkeypatch, tmp_path):
    import firedrake as fd
    cm = calibrate_melt

    def stand_in(name, text):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return str(path)

    msh = stand_in(MESH + ".msh", "$MeshFormat\n")
    obs = stand_in("Melt_Paolo_Davison_Adusumilli_imbie2.csv",
                   "basin,bmr (Gt/yr),bmr uncert (Gt/yr)\n1,10.0,1.0\n2,20.5,2.0\n")
    clim = {v: stand_in(f"meltMIP/OI_Climatology_ismip8km_60m_{v}_extrap.nc", v)
            for v in ("tf", "so")}
    basins = stand_in("basin_numbers_ismip8km_v2.nc", "basins")
    bedmachine = stand_in(
        "NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc", "bed")
    for name, value in (("INV_H5", msh), ("OBS_CSV", obs), ("IMBIE2_NC", basins),
                        ("DATA_ROOT", str(tmp_path))):
        monkeypatch.setattr(cm, name, value)
    monkeypatch.setattr(cm, "_bedmachine_path", lambda: bedmachine)
    monkeypatch.setenv("ISMIP7_SITE", "iu_quartz")
    monkeypatch.setenv("SLURM_JOB_PARTITION", "debug")
    monkeypatch.setenv("SLURM_JOB_ID", "10669722")

    mesh = fd.RectangleMesh(4, 2, 4.0, 2.0)
    g = {"floating": np.arange(16) % 2 == 0}
    record = cm.calibration_record(mesh, g, MESH)
    assert (record["vertices"], record["cells"], record["floating_cells"]) == (15, 16, 8)
    assert record["mesh"] == MESH
    assert record["mesh_file"] == {
        "name": MESH + ".msh", "md5": hashlib.md5(b"$MeshFormat\n").hexdigest(),
        "sha256": file_sha256(msh)}
    assert record["obs_table"] == {
        "name": "Melt_Paolo_Davison_Adusumilli_imbie2.csv",
        "sha256": file_sha256(obs), "total_gtyr": 30.5}
    assert record["inputs_sha256"] == {
        "clim_tf": ["OI_Climatology_ismip8km_60m_tf_extrap.nc", file_sha256(clim["tf"])],
        "clim_so": ["OI_Climatology_ismip8km_60m_so_extrap.nc", file_sha256(clim["so"])],
        "imbie2": ["basin_numbers_ismip8km_v2.nc", file_sha256(basins)]}
    assert record["bedmachine"] == os.path.basename(bedmachine)
    assert (record["site"], record["partition"], record["ranks"], record["job"]) == (
        "iu_quartz", "debug", 1, "10669722")
    assert (record["raster_sample"], record["geometry_space"], record["melt_slope"],
            record["oi_version"]) == ("vertex_front", "dg0", "ant", "30_sep")
    assert record["code"] == "unknown" or len(record["code"]) == 40
    assert str(tmp_path) not in json.dumps(record)
    check_melt_calibration_record(dict(record, K=6.5e-5))

    # started by hand: no site, since a workstation's hostname can name a person
    for key in ("ISMIP7_SITE", "SLURM_JOB_PARTITION", "SLURM_JOB_ID"):
        monkeypatch.delenv(key)
    record = cm.calibration_record(mesh, g, MESH)
    assert (record["site"], record["partition"], record["job"]) == (None, None, None)

    # the release the fit read is the one recorded (select_melt_parameters.py
    # reads ISMIP7_OI_VERSION's)
    nov = {v: stand_in(f"nov/{v}.nc", f"nov {v}") for v in ("tf", "so")}
    monkeypatch.setattr(cm, "_oi_climatology_path",
                        lambda root, v, version: nov[v] if version == "06_nov" else clim[v])
    record = cm.calibration_record(mesh, g, MESH, oi_version="06_nov")
    assert record["oi_version"] == "06_nov"
    assert record["inputs_sha256"]["clim_tf"] == ["tf.nc", file_sha256(nov["tf"])]

    # an input rank 0 cannot read stops the record on every rank
    monkeypatch.setattr(cm, "OBS_CSV", str(tmp_path / "absent.csv"))
    with pytest.raises(FileNotFoundError):
        cm.calibration_record(mesh, g, MESH)
