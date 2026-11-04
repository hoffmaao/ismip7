#!/usr/bin/env python3
r"""One experiment of the ISMIP7 ESM set: the historical or a projection of
a CMIP model on the forcing share, named by the environment.

    ISMIP7_ESM=ACCESS-CM2 ISMIP7_SCENARIO=historical mpiexec -n 32 python scripts/esm_experiment.py
    ISMIP7_ESM=ACCESS-CM2 ISMIP7_SCENARIO=ssp370 mpiexec -n 32 python scripts/esm_experiment.py

The core experiments each have a shim naming their ESM and scenario
(scripts/historical, scripts/projections). An ESM-set member is the same
run under other forcing, so this one driver reads the pair from ISMIP7_ESM
and ISMIP7_SCENARIO and runs experiment.run_core_experiment over the period
the scenario owns (icepack2_tools.experiment_sets). The results are named as
the core runs are, hist_<esm tag> or <scenario>_<esm tag>, and a projection
branches from hist_<esm tag> of the same run tag, so a routine runs the set
exactly as it runs the cores (initialization_routine.sh; the runner
projection.sbatch names this driver for any experiment it has no shim for).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from icepack2_tools.experiment_sets import (  # noqa: E402
    esm_set_counter, experiment_name, experiment_names_esm, scenario_period,
)


def main():
    esm = os.environ.get("ISMIP7_ESM")
    scenario = os.environ.get("ISMIP7_SCENARIO")
    if not esm or not scenario:
        sys.exit("esm_experiment.py needs ISMIP7_ESM (the CMIP model, e.g. ACCESS-CM2) "
                 "and ISMIP7_SCENARIO (historical, ssp370, ssp534-over, ...)")
    t_start, t_end = scenario_period(scenario)
    # A runner names the chain for one model; forcing from another under that
    # name would be a well-formed run of the wrong experiment.
    experiment = os.environ.get("ISMIP7_EXPERIMENT")
    if experiment and not experiment_names_esm(experiment, esm):
        sys.exit(f"ISMIP7_EXPERIMENT={experiment} is not a run of ISMIP7_ESM={esm}")
    counter = esm_set_counter(esm, scenario)
    label = f"ESM set {counter}" if counter else "ESM experiment"
    from experiment import run_core_experiment  # noqa: E402  (loads Firedrake)
    run_core_experiment(
        core=label, title=f"{scenario} with {esm} ({t_start:.0f}-{t_end - 1:.0f})",
        name=experiment_name(esm, scenario), esm=esm, scenario=scenario,
        t_start_default=t_start, t_end_default=t_end,
        restart_from_hist=scenario != "historical",
    )


if __name__ == "__main__":
    main()
