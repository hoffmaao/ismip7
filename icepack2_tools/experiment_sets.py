r"""The ISMIP7 experiment sets beyond the core: the forcing each member runs
under, the period its driver owns, and its counter in the submission names.

The protocol's optional ESM experiments are "centered around SSP370" with
further CMIP models (ISMIP7 cheat sheet, April 2026; ismip.org). The forcing
share carries five further models, each with a historical and an SSP3-7.0
tree downscaled by the energy-balance product dEBM2 (ACCESS-CM2, CanESM5,
GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR), and the SSP5-3.4-overshoot scenario
for the two core models. Each further model's projection branches from its
own historical, as the core experiments do; the overshoot projections branch
from the core historicals (C001, C002), which is why they have no historical
of their own here.

The counters are this repository's: the conventions leave the numbering of
the ESM set to the group, the counter names the experiment in every filename
and directory, and a counter typed by hand for the wrong experiment is a
well-formed submission no checker can catch (the core set's reasoning in
write_ismip7_output.py). So it is fixed here and read from here by the writer
(``set_counter``), by preflight and by the forcing mirror.
"""

#: The CMIP models of the ESM set that have no core experiment.
ESM_SET_MODELS = ("ACCESS-CM2", "CanESM5", "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR")

#: ``(counter, esm, scenario)`` of every ESM-set member, in counter order.
ESM_SET = (
    ("E001", "ACCESS-CM2", "historical"),
    ("E002", "ACCESS-CM2", "ssp370"),
    ("E003", "CanESM5", "historical"),
    ("E004", "CanESM5", "ssp370"),
    ("E005", "GFDL-ESM4", "historical"),
    ("E006", "GFDL-ESM4", "ssp370"),
    ("E007", "IPSL-CM6A-LR", "historical"),
    ("E008", "IPSL-CM6A-LR", "ssp370"),
    ("E009", "MPI-ESM1-2-HR", "historical"),
    ("E010", "MPI-ESM1-2-HR", "ssp370"),
    ("E011", "CESM2-WACCM", "ssp534-over"),
    ("E012", "MRI-ESM2-0", "ssp534-over"),
)

#: The scenarios a forward can be forced with, and the period each owns, as
#: the core shims declare it: ``t=Y.0`` is 1 January of Y, so a historical
#: covering 2003 to 2014 ends at 2015, SSP3-7.0 (to 2100) at 2101 and a
#: scenario to 2300 at 2301.
SCENARIO_PERIODS = {
    "historical": (2003.0, 2015.0),
    "ssp126": (2015.0, 2301.0),
    "ssp370": (2015.0, 2101.0),
    "ssp534-over": (2015.0, 2301.0),
    "ssp585": (2015.0, 2301.0),
}

# The experiment suffixes the core runners use for the two core models, which
# predate the rule below (projection.sbatch, initialization_routine.sh).
_CORE_KEYS = {"CESM2-WACCM": "cesm_waccm", "MRI-ESM2-0": "mri_esm2"}


def esm_set_counter(esm, scenario):
    r"""The ESM set's counter of ``(esm, scenario)``, e.g. ``E002``, or None
    when the pair is not a member."""
    for counter, e, s in ESM_SET:
        if (e, s) == (esm, scenario):
            return counter
    return None


def esm_set_members():
    r"""``[(esm, scenario)]`` in counter order."""
    return [(e, s) for _, e, s in ESM_SET]


def scenario_period(scenario):
    r"""``(t_start, t_end)`` the drivers default to for ``scenario``."""
    try:
        return SCENARIO_PERIODS[scenario]
    except KeyError:
        raise ValueError(f"no experiment period for scenario {scenario!r}; one of "
                         f"{', '.join(SCENARIO_PERIODS)}") from None


def esm_tag(esm):
    r"""The ESM's element of a results name: ``MRI-ESM2-0`` is ``mri_esm2_0``,
    as experiment.py builds it."""
    return esm.lower().replace("-", "_")


def experiment_name(esm, scenario):
    r"""The results stem of a run (before the run tag and the resolution), as
    the core shims name theirs: ``hist_<tag>`` or ``<scenario>_<tag>``."""
    head = "hist" if scenario == "historical" else scenario
    return f"{head}_{esm_tag(esm)}"


def experiment_key(esm):
    r"""The ESM's suffix in a runner's experiment name (``ISMIP7_EXPERIMENT``,
    a routine's stage name): the core pair keep the suffixes the runners were
    written with, every other model is named by its tag."""
    return _CORE_KEYS.get(esm, esm_tag(esm))


def experiment_names_esm(experiment, esm):
    r"""Whether the runner's experiment name ``experiment`` (``hist_<key>`` or
    ``<scenario>_<key>``) is a run of ``esm``: the one driver of the ESM set
    takes its model from ``ISMIP7_ESM`` and refuses a chain named for
    another, so ACCESS-CM2 forcing never runs under a CanESM5 name."""
    return experiment.endswith("_" + experiment_key(esm))


def experiment_scenario(experiment):
    r"""The scenario of a runner's experiment name: ``historical`` for
    ``hist_<key>``, else the part before the first underscore."""
    if experiment.startswith("hist_"):
        return "historical"
    return experiment.split("_", 1)[0]
