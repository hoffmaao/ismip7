#!/usr/bin/env python3
r"""The relaxation year of the relaxed initial state (icepack2_tools.relaxation).

Cold-starts from the MAP that ``ISMIP7_INVERSION`` names one year before the
geometry year (``ISMIP7_RELAX_START``, default 2014), with that year of the
Smith et al. (2020) mean dH/dt undone on grounded ice (issue #117), and runs
to 2015.0 at ``ISMIP7_RELAX_DT`` (default half of ``ISMIP7_DT``):

    forcing         OCX's own for the year, as ISMIP7_OCX_FORCING selects it
                    (ISMIP7_RELAX_FORCING=none: no SMB or melt, a smoke test
                    that needs no forcing data)
    apparent MB     off, so the geometry answers the model's own imbalance
    front           pinned (ISMIP7_FIXED_FRONT=1, ISMIP7_CALVING=none), so the
                    re-inversion keeps the MAP's extent

The final checkpoint, ``results/relax_<MAP stem>_<lc>_final.h5``, seeds the
re-inversion (``ISMIP7_WARM_START`` and ``ISMIP7_MESH=checkpoint``). It records
the MAP, its sha256, the year, the step and the forcing, and carries the
MAP's objective settings, so the re-inversion minimises the MAP's objective
on the relaxed geometry. Everything else the driver needs is a forward's.

Usage:
    submit.sh projection ISMIP7_EXPERIMENT=relax ISMIP7_INVERSION=<MAP> ISMIP7_MESH=checkpoint
    mpiexec -n 24 python scripts/relaxation/run.py
"""
import os
import sys

_SCRIPTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from mpi4py import MPI  # noqa: E402

from simulation import (setup_model, run_simulation, auto_resume_checkpoint,  # noqa: E402
                        auto_resume, PETSc)
from projections.ocx import (protocol_forcing, ocx_feedback,  # noqa: E402
                             ocx_forcing_callback)
from icepack2_tools.forcing import feedback_mode, reject_collapse_mask  # noqa: E402
from icepack2_tools.handoff import OBJECTIVE_KEYS, OBJECTIVE_RECORD_KEYS  # noqa: E402
from icepack2_tools.relaxation import (  # noqa: E402
    END_STATE_ATTR, describe_relaxation, relax_backdate_years, relax_dt,
    relax_environment_problems, relax_experiment_name, relax_forcing,
    relax_window, relaxation_attrs,
)
from icepack2_tools.runconfig import (  # noqa: E402
    apparent_mb_mode, calving_law, deltat_per_basin_npz, file_sha256,
    fixed_front, ismip7_output, ocx_forcing,
)

OUTPUT_INTERVAL = int(os.environ.get("ISMIP7_OUTPUT_INTERVAL", "10"))


def _root_attrs(path, keys):
    r"""The root attributes ``keys`` of the HDF5 file ``path`` that it has,
    read on rank 0 and shared."""
    comm = MPI.COMM_WORLD
    attrs = None
    if comm.rank == 0:
        import h5py
        with h5py.File(path, "r") as h:
            attrs = {k: h["/"].attrs[k] for k in keys if k in h["/"].attrs}
    return comm.bcast(attrs, root=0)


def main():
    comm = MPI.COMM_WORLD
    source = os.environ.get("ISMIP7_INVERSION", "").strip()
    if not source:
        raise RuntimeError("ISMIP7_INVERSION must name the MAP to relax")
    if not os.path.exists(source):
        raise FileNotFoundError(f"ISMIP7_INVERSION={source} does not exist")
    t_start, t_end = relax_window()
    step = relax_dt(t_start, t_end)
    backdate = relax_backdate_years(t_start)
    forcing = relax_forcing()
    problems = relax_environment_problems(
        apparent_mb=apparent_mb_mode(), output=ismip7_output(),
        calving=calving_law(), fixed_front=fixed_front())
    if problems:
        raise RuntimeError("The relaxation year cannot run like this: "
                           + " ".join(problems))
    reject_collapse_mask("the relaxation year")

    source_attrs = _root_attrs(
        source, OBJECTIVE_KEYS + OBJECTIVE_RECORD_KEYS
        + ("geometry_source_method", END_STATE_ATTR))
    sha = comm.bcast(file_sha256(source) if comm.rank == 0 else None, root=0)
    forcing_name = (f"ocx {ocx_forcing()}" if forcing == "ocx" else "none")
    # Built before the setup: a relaxed MAP is refused here, not after it.
    final_attrs = relaxation_attrs(
        source_map=source, source_sha256=sha, source_attrs=source_attrs,
        t_start=t_start, t_end=t_end, dt=step, backdate_years=backdate,
        forcing=forcing_name)

    experiment_name = relax_experiment_name(source, os.environ.get("ISMIP7_RUN_TAG", ""))
    # A link the wall clock stopped saves its final state with the record
    # above, so a resume continues the same year of the same MAP; a state
    # from another MAP under the same name is refused.
    restart = os.environ.get("ISMIP7_RESTART")
    if restart is None and auto_resume():
        restart = auto_resume_checkpoint(experiment_name)
        PETSc.Sys.Print(f"Auto-resume: {restart}" if restart
                        else "Auto-resume: no prior checkpoint")
    if restart:
        _had = _root_attrs(restart, ("relax_source_sha256",)).get("relax_source_sha256")
        _had = _had.decode() if isinstance(_had, bytes) else _had
        if _had != sha:
            raise RuntimeError(
                f"{restart} is not a relaxation of {source} (sha256 {_had} against "
                f"{sha}); move it aside to relax this MAP from the start")

    readers = None
    feedback = None
    dT_npz = deltat_per_basin_npz()
    if forcing == "ocx":
        readers = (protocol_forcing(t_start, t_end)
                   if ocx_forcing() == "protocol" else None)
        feedback = ocx_feedback(readers, t_start, t_end)

    ctx = setup_model(
        restart_from=restart,
        backdate_years=0.0 if restart else backdate,
        smb_feedback=feedback_mode(feedback) if forcing == "ocx" else None)
    if not restart and not ctx["map_geometry_taken"]:
        raise RuntimeError(
            f"the relaxation year runs on the mesh and geometry of {source}, "
            f"and this run built its own geometry on {ctx['mesh_basename']}; "
            f"set ISMIP7_MESH=checkpoint")

    callback, ocean = None, None
    if forcing == "ocx":
        callback, lines, ocean = ocx_forcing_callback(ctx, readers, feedback, dT_npz)
        for line in lines:
            PETSc.Sys.Print(f"  {line}")
    else:
        PETSc.Sys.Print("  Forcing: none (ISMIP7_RELAX_FORCING=none): no SMB, no melt")
    ctx["final_attrs"] = final_attrs

    PETSc.Sys.Print("\nRelaxation year of the relaxed initial state")
    PETSc.Sys.Print(f"  {describe_relaxation(final_attrs)}")
    PETSc.Sys.Print(f"  source sha256 {sha}")
    PETSc.Sys.Print(f"  Period: {t_start}-{t_end}, dt={step} "
                    f"({int(round((t_end - t_start) / step))} steps), "
                    f"{backdate:g} yr of Smith dH/dt undone at the start")

    run_simulation(
        ctx,
        experiment_name=experiment_name,
        t_start=t_start,
        t_end=t_end,
        dt=step,
        output_interval=OUTPUT_INTERVAL,
        forcing_callback=callback,
    )
    if ocean is not None:
        ocean.close()


if __name__ == "__main__":
    main()
