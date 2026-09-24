r"""An unattended run resumes only a checkpoint written under its own law.

The drivers auto-resume from the newest checkpoint of their experiment name,
which carries only the run tag. A law-driven run without a tag of its own
would otherwise pick up a stock run's state (or a stock run a law run's) and
then overwrite that run's files. ``simulation.auto_resume_checkpoint`` is
what all three drivers now ask.
"""
import importlib.util
import os

import pytest
from firedrake import CheckpointFile, UnitSquareMesh

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "antarctica", "scripts")
NAME = "ctrl2015_cesm2_waccm"


@pytest.fixture
def simulation(tmp_path, monkeypatch):
    import sys
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    for knob in ("ISMIP7_CALVING", "ISMIP7_CALVING_PARAMS", "ISMIP7_CALVING_MODULE"):
        monkeypatch.delenv(knob, raising=False)
    spec = importlib.util.spec_from_file_location(
        "simulation", os.path.join(SCRIPTS, "simulation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS_DIR", str(tmp_path))
    return mod


def _checkpoint(sim, law=None, t_yr=2030.0):
    path = os.path.join(sim.RESULTS_DIR, f"{NAME}_1000_final.h5")
    with CheckpointFile(path, "w") as chk:
        chk.save_mesh(UnitSquareMesh(2, 2))
        chk.set_attr("/", "t_yr", t_yr)
        if law is not None:
            chk.set_attr("/", "calving_law", law)
    return path


def _law(monkeypatch, name, params=""):
    monkeypatch.setenv("ISMIP7_CALVING", name)
    if params:
        monkeypatch.setenv("ISMIP7_CALVING_PARAMS", params)
    from icepack2_tools.runconfig import calving_law_object
    return calving_law_object().describe()


def test_nothing_to_resume_is_none(simulation):
    assert simulation.auto_resume_checkpoint(NAME, 1000) is None


def test_a_stock_run_resumes_its_own_checkpoint(simulation):
    path = _checkpoint(simulation)
    assert simulation.auto_resume_checkpoint(NAME, 1000) == path


def test_a_law_run_resumes_its_own_checkpoint(simulation, monkeypatch):
    law = _law(monkeypatch, "vonmises", "sigma_max_fl=0.2")
    path = _checkpoint(simulation, law)
    assert simulation.auto_resume_checkpoint(NAME, 1000) == path


def test_a_law_run_does_not_resume_a_stock_run(simulation, monkeypatch):
    _checkpoint(simulation)
    _law(monkeypatch, "vonmises")
    with pytest.raises(RuntimeError, match="calving law none, but this run is configured with vonmises"):
        simulation.auto_resume_checkpoint(NAME, 1000)


def test_a_stock_run_does_not_resume_a_law_run(simulation, monkeypatch):
    law = _law(monkeypatch, "thickness", "hc=150")
    _checkpoint(simulation, law)
    monkeypatch.delenv("ISMIP7_CALVING")
    monkeypatch.delenv("ISMIP7_CALVING_PARAMS")
    with pytest.raises(RuntimeError, match="configured with none"):
        simulation.auto_resume_checkpoint(NAME, 1000)


def test_a_changed_parameter_is_a_different_run(simulation, monkeypatch):
    _checkpoint(simulation, _law(monkeypatch, "vonmises", "sigma_max_fl=0.2"))
    monkeypatch.setenv("ISMIP7_CALVING_PARAMS", "sigma_max_fl=0.25")
    with pytest.raises(RuntimeError, match="ISMIP7_RUN_TAG"):
        simulation.auto_resume_checkpoint(NAME, 1000)
