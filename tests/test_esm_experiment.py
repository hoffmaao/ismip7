r"""The ESM set's one driver (``antarctica/scripts/esm_experiment.py``)
checks its environment before it loads Firedrake: a chain named for another
scenario is refused, and a further model's aSMB pool is historical + ssp370,
the one projection it carries, while the core models keep the shared
default. The run itself is replaced by a recorder."""

import importlib.util
import os
import sys
import types

import pytest

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "antarctica", "scripts", "esm_experiment.py")


@pytest.fixture
def driver(monkeypatch):
    calls, printed = [], []

    def run_core_experiment(**kwargs):
        calls.append((kwargs, os.environ.get("ISMIP7_CLIM_SCENARIO")))

    monkeypatch.setitem(sys.modules, "experiment",
                        types.SimpleNamespace(run_core_experiment=run_core_experiment))
    petsc = types.SimpleNamespace(Sys=types.SimpleNamespace(Print=printed.append))
    monkeypatch.setitem(sys.modules, "firedrake", types.ModuleType("firedrake"))
    monkeypatch.setitem(sys.modules, "firedrake.petsc", types.SimpleNamespace(PETSc=petsc))
    for v in ("ISMIP7_ESM", "ISMIP7_SCENARIO", "ISMIP7_EXPERIMENT", "ISMIP7_CLIM_SCENARIO"):
        monkeypatch.delenv(v, raising=False)
    spec = importlib.util.spec_from_file_location("esm_experiment_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def run(**env):
        for k, v in env.items():
            monkeypatch.setenv(k, v)
        module.main()
        return calls, printed
    return run


def test_a_further_model_pools_with_ssp370_by_default(driver):
    calls, printed = driver(ISMIP7_ESM="ACCESS-CM2", ISMIP7_SCENARIO="ssp370",
                            ISMIP7_EXPERIMENT="ssp370_access_cm2")
    (kwargs, pool), = calls
    assert pool == "ssp370"
    assert kwargs["name"] == "ssp370_access_cm2" and kwargs["core"] == "ESM set E002"
    assert printed and "ISMIP7_CLIM_SCENARIO=ssp370" in printed[0]


def test_a_further_model_refuses_a_pool_it_does_not_carry(driver):
    with pytest.raises(SystemExit, match="ISMIP7_CLIM_SCENARIO=ssp126"):
        driver(ISMIP7_ESM="CanESM5", ISMIP7_SCENARIO="historical",
               ISMIP7_CLIM_SCENARIO="ssp126")


def test_a_core_model_keeps_the_shared_default_pool(driver):
    calls, printed = driver(ISMIP7_ESM="CESM2-WACCM", ISMIP7_SCENARIO="ssp534-over",
                            ISMIP7_EXPERIMENT="ssp534-over_cesm_waccm")
    (kwargs, pool), = calls
    assert pool is None and printed == []
    assert kwargs["core"] == "ESM set E011"


def test_a_chain_named_for_another_scenario_is_refused(driver):
    with pytest.raises(SystemExit, match="not a run of ISMIP7_SCENARIO=historical"):
        driver(ISMIP7_ESM="ACCESS-CM2", ISMIP7_SCENARIO="historical",
               ISMIP7_EXPERIMENT="ssp370_access_cm2")
