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
    ESM_SET_MODELS, esm_set_counter, experiment_name, experiment_names_esm,
    experiment_scenario, scenario_period,
)

#: The one projection the further models carry on the share, and so the one
#: scenario their aSMB pool (``ISMIP7_CLIM_SCENARIO``) can pool with.
FURTHER_MODEL_POOL = "ssp370"


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
    if experiment and experiment_scenario(experiment) != scenario:
        sys.exit(f"ISMIP7_EXPERIMENT={experiment} is not a run of ISMIP7_SCENARIO={scenario}")
    # A further model has no ssp126 tree, so the shared default pool would
    # re-reference its aSMB over the historical years alone; the core models
    # keep the shared default, so E011 and E012 pool as the cores do.
    pool_note = None
    if esm in ESM_SET_MODELS:
        pool = os.environ.get("ISMIP7_CLIM_SCENARIO")
        if pool is None:
            os.environ["ISMIP7_CLIM_SCENARIO"] = FURTHER_MODEL_POOL
            pool_note = f"ISMIP7_CLIM_SCENARIO={FURTHER_MODEL_POOL}: the one projection {esm} carries"
        elif pool != FURTHER_MODEL_POOL:
            sys.exit(f"ISMIP7_CLIM_SCENARIO={pool} names a scenario {esm} does not carry; "
                     f"its aSMB pool is historical + {FURTHER_MODEL_POOL}")
    counter = esm_set_counter(esm, scenario)
    label = f"ESM set {counter}" if counter else "ESM experiment"
    from experiment import run_core_experiment  # noqa: E402  (loads Firedrake)
    if pool_note:
        from firedrake.petsc import PETSc
        PETSc.Sys.Print(pool_note)
    run_core_experiment(
        core=label, title=f"{scenario} with {esm} ({t_start:.0f}-{t_end - 1:.0f})",
        name=experiment_name(esm, scenario), esm=esm, scenario=scenario,
        t_start_default=t_start, t_end_default=t_end,
        restart_from_hist=scenario != "historical",
    )


if __name__ == "__main__":
    main()
