r"""The ESM set's table (``icepack2_tools.experiment_sets``): the counters
that name its members in every submission filename, the periods its one
driver runs, and the experiment names a runner and the routine build from
a model's name, which have to agree with the core runners' for the two core
models and with ``experiment.py``'s results names for every model."""

import pytest

from icepack2_tools.experiment_sets import (
    ESM_SET, ESM_SET_MODELS, clim_scenario_for, esm_set_counter, esm_set_members, esm_tag,
    experiment_key, experiment_name, experiment_names_esm, experiment_scenario,
    scenario_period,
)


def test_the_counters_are_consecutive_and_name_each_member_once():
    counters = [c for c, _, _ in ESM_SET]
    assert counters == [f"E{n:03d}" for n in range(1, len(ESM_SET) + 1)]
    assert len(set(esm_set_members())) == len(ESM_SET)


def test_every_further_model_has_a_historical_and_its_projection():
    r"""The protocol runs one historical per CMIP model and the ESM
    experiments centre on SSP3-7.0, so each further model is two members;
    the overshoot projections of the core models branch from the core
    historicals and are one member each."""
    members = set(esm_set_members())
    for esm in ESM_SET_MODELS:
        assert (esm, "historical") in members and (esm, "ssp370") in members
    assert ("CESM2-WACCM", "ssp534-over") in members and ("MRI-ESM2-0", "ssp534-over") in members
    assert ("CESM2-WACCM", "historical") not in members


def test_a_counter_follows_from_the_forcing():
    assert esm_set_counter("ACCESS-CM2", "historical") == "E001"
    assert esm_set_counter("ACCESS-CM2", "ssp370") == "E002"
    assert esm_set_counter("MPI-ESM1-2-HR", "ssp370") == "E010"
    assert esm_set_counter("MRI-ESM2-0", "ssp534-over") == "E012"
    assert esm_set_counter("CESM2-WACCM", "ssp370") is None, "a core experiment is not in the set"


def test_the_period_is_the_scenario_s_own():
    r"""``t=Y.0`` is 1 January of Y, as the core shims declare their ends."""
    assert scenario_period("historical") == (2003.0, 2015.0)
    assert scenario_period("ssp370") == (2015.0, 2101.0)
    assert scenario_period("ssp534-over") == (2015.0, 2301.0)
    with pytest.raises(ValueError, match="no experiment period"):
        scenario_period("ssp245")


def test_results_names_follow_experiment_py():
    r"""``experiment.py`` names a run ``<name>_<tag>_<lc>`` with the ESM's
    tag ``lower().replace('-', '_')``, and ``historical_endpoint`` looks a
    projection's start up under ``hist_<tag>``; the ESM set's names are the
    same construction, so the core shims' names are reproduced exactly."""
    assert esm_tag("MRI-ESM2-0") == "mri_esm2_0"
    assert experiment_name("MRI-ESM2-0", "historical") == "hist_mri_esm2_0"
    assert experiment_name("CESM2-WACCM", "ssp370") == "ssp370_cesm2_waccm"
    assert experiment_name("ACCESS-CM2", "ssp370") == "ssp370_access_cm2"
    assert experiment_name("IPSL-CM6A-LR", "historical") == "hist_ipsl_cm6a_lr"


def test_runner_names_keep_the_core_suffixes_and_tag_the_rest():
    r"""projection.sbatch and the routine name the core pair ``cesm_waccm``
    and ``mri_esm2`` (older than this table); a further model is named by
    its tag, and a name is a run of exactly one model."""
    assert experiment_key("CESM2-WACCM") == "cesm_waccm"
    assert experiment_key("MRI-ESM2-0") == "mri_esm2"
    assert experiment_key("MPI-ESM1-2-HR") == "mpi_esm1_2_hr"
    assert experiment_names_esm("ssp370_access_cm2", "ACCESS-CM2")
    assert experiment_names_esm("hist_mri_esm2", "MRI-ESM2-0")
    assert not experiment_names_esm("ssp370_access_cm2", "CanESM5")
    assert not experiment_names_esm("ssp370_canesm5", "ACCESS-CM2")
    assert experiment_scenario("hist_access_cm2") == "historical"
    assert experiment_scenario("ssp534-over_cesm_waccm") == "ssp534-over"
    assert experiment_scenario("ssp370_mpi_esm1_2_hr") == "ssp370"


def test_a_further_model_pools_with_its_one_projection():
    from icepack2_tools.climatology import CLIM_SCENARIO_DEFAULT
    assert clim_scenario_for("CanESM5", None) == "ssp370"
    assert clim_scenario_for("CanESM5", "ssp370") == "ssp370"
    with pytest.raises(ValueError, match="does not carry"):
        clim_scenario_for("CanESM5", "ssp126")
    assert clim_scenario_for("CESM2-WACCM", None) == CLIM_SCENARIO_DEFAULT
    assert clim_scenario_for("MRI-ESM2-0", "ssp370") == "ssp370"
