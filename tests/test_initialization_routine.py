r"""The staged initialization routine's state machine.

``antarctica/scripts/initialization_routine.sh`` turns the document's stages
into submissions of the batch runners, each under the environment it recorded
at ``init`` plus the stage's own knobs, and hands itself on from the link that
finishes a stage. What it must get right is the order, the inputs (which MAP
or state each stage starts from, under which name), the knobs that may not
leak between stages, and idempotence: a call that finds nothing new to start
starts nothing. Those decisions run here against a stub ``sbatch`` that
records what it was given, with the stage products faked on disk.
"""

import os
import subprocess
import sys
from pathlib import Path

import h5py
import pytest

REPO = Path(__file__).resolve().parent.parent

SBATCH_STUB = (
    "#!/bin/bash\n"
    'printf "ARGV: %s\\n" "$*" >> "$SBATCH_CALLS"\n'
    'env | grep -E "^(ISMIP7_|SLURM_)" | sort | sed "s/^/ENV: /" >> "$SBATCH_CALLS"\n'
    'echo "END" >> "$SBATCH_CALLS"\n'
    'echo "$((RANDOM + 1000))"\n'
)


@pytest.fixture
def sandbox(tmp_path):
    r"""A checkout whose scripts are the repository's and whose results
    directory is the sandbox's own."""
    ant = tmp_path / "antarctica"
    ant.mkdir()
    (ant / "scripts").symlink_to(REPO / "antarctica" / "scripts")
    (tmp_path / "icepack2_tools").symlink_to(REPO / "icepack2_tools")
    (ant / "results").mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "sbatch").write_text(SBATCH_STUB)
    (bin_dir / "sbatch").chmod(0o755)
    (bin_dir / "python").write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    (bin_dir / "python").chmod(0o755)
    (tmp_path / "out").mkdir()
    return tmp_path


def base_env(sandbox):
    calls = sandbox / "sbatch_calls.txt"
    return {
        "PATH": f"{sandbox / 'bin'}:{os.environ['PATH']}",
        "HOME": str(sandbox),
        "SBATCH_CALLS": str(calls),
        "USER": "tester",
        # the objective the routine records
        "ISMIP7_LC": "32000",
        "ISMIP7_MESH": "/meshes/antarctica_320000_32000.msh",
        "ISMIP7_MISFIT_SCALE": "14000",
        "ISMIP7_DHDT_WEIGHT": "1",
        "ISMIP7_GRAD_PRECOND": "mass_consistent",
        "PYTHONPATH": str(sandbox),
        # a stage knob left in the caller's shell must not be recorded
        "ISMIP7_INVERT": "both",
        "ISMIP7_APPARENT_MB": "0",
        "ISMIP7_ROUTINE_LAUNCH": "slurm",
        "ISMIP7_ROUTINE_SBATCH_INV": "-p scavenge --time=01:00:00 --mem=240G",
        "ISMIP7_ROUTINE_SBATCH_FWD": "-p scavenge --time=01:00:00 --mem=180G",
        "ISMIP7_ROUTINE_ITERS": "40",
        "ISMIP7_ROUTINE_RELAX_YR": "20",
        "ISMIP7_ROUTINE_WARM_START": str(sandbox / "warm.h5"),
    }


def routine(sandbox, *args, env=None):
    merged = base_env(sandbox)
    merged.update(env or {})
    script = sandbox / "antarctica" / "scripts" / "initialization_routine.sh"
    proc = subprocess.run(
        ["bash", str(script), str(sandbox / "out"), *args], env=merged,
        cwd=str(sandbox), capture_output=True, text=True, timeout=120)
    return proc.returncode, proc.stdout + proc.stderr


def submits(sandbox):
    r"""One dict per recorded sbatch call: argv and the ISMIP7 environment,
    then the record is cleared."""
    calls = sandbox / "sbatch_calls.txt"
    if not calls.exists():
        return []
    out = []
    for block in calls.read_text().split("END\n"):
        if not block.strip():
            continue
        lines = block.strip().splitlines()
        env = dict(line[5:].split("=", 1) for line in lines[1:] if line.startswith("ENV: "))
        out.append({"argv": lines[0][6:], "env": env})
    calls.unlink()
    return out


def finish_inversion(path, its=40, **attrs):
    with h5py.File(path, "w") as h:
        h["/"].attrs["objective_iteration"] = its
        h["/"].attrs["objective_total"] = 1.0e6
        h["/"].attrs["misfit_scale"] = 14000.0
        h["/"].attrs["fluidity_control"] = "floating"
        for key, value in attrs.items():
            h["/"].attrs[key] = value
    Path(str(path) + ".done").touch()


def finish_forward(path, t_yr):
    with h5py.File(path, "w") as h:
        h["/"].attrs["t_yr"] = float(t_yr)
        h["/"].attrs["dt_yr"] = 0.05
        h["/"].attrs["stalled"] = 0
        h["/"].attrs["fssa_reference"] = "none"


def by_name(subs, name):
    found = [s for s in subs if f"-J {name} " in s["argv"] or s["argv"].endswith(f"-J {name}")]
    assert len(found) == 1, f"{name}: {[s['argv'] for s in subs]}"
    return found[0]


def test_the_routine_runs_itself_stage_by_stage(sandbox):
    out = sandbox / "out"
    results = sandbox / "antarctica" / "results"
    (sandbox / "warm.h5").write_text("warm\n")

    rc, log = routine(sandbox, "init")
    assert rc == 0, log
    env_text = (out / "routine.env").read_text()
    assert 'export ISMIP7_MISFIT_SCALE="14000"' in env_text
    assert "ISMIP7_INVERT" not in env_text, "stage knobs are never recorded"
    assert "ISMIP7_APPARENT_MB" not in env_text
    assert 'ISMIP7_ROUTINE_TAG=out' in env_text

    # init starts the first stage and nothing else: everything after it waits
    subs = submits(sandbox)
    assert [s["argv"].split(" -J ")[1].split()[0] for s in subs] == ["rt_out_shelf1"]
    s1 = subs[0]["env"]
    assert s1["ISMIP7_INVERT"] == "phi"
    assert s1["ISMIP7_PHI_BOUNDS"] == "-2.3,2.3"
    assert s1["ISMIP7_FLUIDITY_CONTROL"] == "floating"
    assert s1["ISMIP7_MAXITER"] == "40" and s1["ISMIP7_CHAIN_ITER_BUDGET"] == "40"
    assert s1["ISMIP7_MAP_OUT"] == str(out / "map_shelf1.h5")
    assert s1["ISMIP7_WARM_START"] == str(sandbox / "warm.h5")
    assert s1["ISMIP7_WARM_START_STRICT"] == "0"
    assert s1["ISMIP7_CHAIN_THEN_CMD"].endswith(f"initialization_routine.sh {out} next")
    assert "ISMIP7_DHDT_REWEIGHT" not in s1
    assert "SLURM_JOB_ID" not in s1
    assert "inversion.sbatch" in subs[0]["argv"]
    assert "--time=01:00:00 --mem=240G" in subs[0]["argv"]

    # a second call finds the stage started and starts nothing
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    assert submits(sandbox) == []

    # the friction stage warm-starts from the shelf MAP and re-weights the fit
    finish_inversion(out / "map_shelf1.h5")
    rc, log = routine(sandbox, "next", env={"ISMIP7_INVERT": "phi", "SLURM_JOB_ID": "77"})
    assert rc == 0, log
    (s2,) = submits(sandbox)
    assert "-J rt_out_friction1" in s2["argv"]
    assert s2["env"]["ISMIP7_INVERT"] == "theta"
    assert s2["env"]["ISMIP7_DHDT_REWEIGHT"] == "24"
    assert "ISMIP7_PHI_BOUNDS" not in s2["env"]
    assert s2["env"]["ISMIP7_WARM_START"] == str(out / "map_shelf1.h5")
    assert s2["env"]["ISMIP7_DHDT_WEIGHT"] == "1", "the recorded objective reaches every stage"
    assert "SLURM_JOB_ID" not in s2["env"], "the calling job's scheduler variables stay behind"

    finish_inversion(out / "map_friction1.h5")
    routine(sandbox, "next")
    (s3,) = submits(sandbox)
    assert "-J rt_out_shelf2" in s3["argv"] and s3["env"]["ISMIP7_MAXITER"] == "20"
    finish_inversion(out / "map_shelf2.h5")
    routine(sandbox, "next")
    (s4,) = submits(sandbox)
    assert "-J rt_out_friction2" in s4["argv"] and s4["env"]["ISMIP7_CHAIN_ITER_BUDGET"] == "20"

    # the relaxation: the control driver, apparent MB off, front pinned
    finish_inversion(out / "map_friction2.h5", its=20, misfit_norm="sigma")
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    (relax,) = submits(sandbox)
    assert "-J rt_out_relax" in relax["argv"] and "projection.sbatch" in relax["argv"]
    e = relax["env"]
    assert e["ISMIP7_EXPERIMENT"] == "control" and e["ISMIP7_ESM"] == "CESM2-WACCM"
    assert e["ISMIP7_INVERSION"] == str(out / "map_friction2.h5")
    assert e["ISMIP7_RUN_TAG"] == "out_relax" and e["ISMIP7_T_END"] == "2035"
    assert e["ISMIP7_APPARENT_MB"] == "0" and e["ISMIP7_OUTPUT"] == "0"
    assert e["ISMIP7_FIXED_FRONT"] == "1" and e["ISMIP7_CALVING"] == "none"
    assert "ISMIP7_INVERT" not in e and "ISMIP7_MAP_OUT" not in e
    assert "--mem=180G" in relax["argv"]

    # short of its end year the relaxation is still in progress
    final = results / "ctrl2015_cesm2_waccm_out_relax_32000_final.h5"
    finish_forward(final, 2030)
    routine(sandbox, "next")
    assert submits(sandbox) == []
    rc, log = routine(sandbox, "status")
    assert rc == 0 and "relax                in progress (t=2030.0 of 2035)" in log, log

    # at its end the state is stamped and the re-inversion starts from it,
    # under the name the relaxation module requires
    finish_forward(final, 2035)
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    state = out / "relaxed_state.h5"
    assert state.exists()
    with h5py.File(state, "r") as h:
        a = h["/"].attrs
        assert int(a["relaxation_end_state"]) == 1
        assert float(a["relax_t_start"]) == 2015.0 and float(a["relax_t_end"]) == 2035.0
        assert os.path.basename(a["relax_source_map"]) == "map_friction2.h5"
        assert float(a["misfit_scale"]) == 14000.0, "the MAP's objective travels with the state"
        assert float(a["relax_dt"]) == 0.05
    (r1,) = submits(sandbox)
    assert "-J rt_out_r_shelf1" in r1["argv"]
    assert r1["env"]["ISMIP7_WARM_START"] == str(state)
    assert r1["env"]["ISMIP7_MAP_OUT"] == str(out / "map_friction2_relax2015.h5")
    assert r1["env"]["ISMIP7_INVERT"] == "phi"

    finish_inversion(out / "map_friction2_relax2015.h5")
    routine(sandbox, "next")
    (r2,) = submits(sandbox)
    assert "-J rt_out_r_friction1" in r2["argv"]
    assert r2["env"]["ISMIP7_WARM_START"] == str(out / "map_friction2_relax2015.h5")
    finish_inversion(out / "map_r_friction1.h5")
    routine(sandbox, "next")
    (r3,) = submits(sandbox)
    assert "-J rt_out_r_shelf2" in r3["argv"]
    finish_inversion(out / "map_r_shelf2.h5")
    routine(sandbox, "next")
    (r4,) = submits(sandbox)
    assert "-J rt_out_r_friction2" in r4["argv"]

    # the last MAP is the forward MAP; the historicals (one per ESM) and OCX
    # start from it with the apparent MB off and the submission's front
    finish_inversion(out / "map_r_friction2.h5", its=20)
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    assert os.readlink(out / "map_forward.h5") == "map_r_friction2.h5"
    subs = submits(sandbox)
    names = sorted(s["argv"].split(" -J ")[1].split()[0] for s in subs)
    assert names == ["rt_out_hist_cesm_waccm", "rt_out_hist_mri_esm2", "rt_out_ocx"]
    h = by_name(subs, "rt_out_hist_mri_esm2")["env"]
    assert h["ISMIP7_EXPERIMENT"] == "hist_mri_esm2" and h["ISMIP7_ESM"] == "MRI-ESM2-0"
    assert h["ISMIP7_INVERSION"] == str(out / "map_forward.h5")
    assert h["ISMIP7_RUN_TAG"] == "out" and h["ISMIP7_T_START"] == "2002"
    assert h["ISMIP7_APPARENT_MB"] == "0" and h["ISMIP7_OUTPUT"] == "1"
    assert h["ISMIP7_CALVING"] == "hfb" and h["ISMIP7_FRONT_ADVANCE"] == "none"
    assert h["ISMIP7_FRACTURE"] == "none"
    assert h["ISMIP7_GEOMETRY_BACKDATE"] == "12", "the forced 2002 year starts from the 2003 geometry"
    assert "ISMIP7_FIXED_FRONT" not in h and "ISMIP7_T_END" not in h
    assert "ISMIP7_INVERT" not in h and "ISMIP7_FLUIDITY_CONTROL" not in h

    # one historical at its end releases its own control and projections only
    finish_forward(results / "hist_cesm2_waccm_out_32000_final.h5", 2015)
    finish_forward(results / "hist_mri_esm2_0_out_32000_final.h5", 2010)
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    subs = submits(sandbox)
    names = sorted(s["argv"].split(" -J ")[1].split()[0] for s in subs)
    assert names == ["rt_out_control_cesm_waccm", "rt_out_ssp126_cesm_waccm",
                     "rt_out_ssp370_cesm_waccm", "rt_out_ssp585_cesm_waccm"]
    c = by_name(subs, "rt_out_control_cesm_waccm")["env"]
    assert c["ISMIP7_EXPERIMENT"] == "control" and c["ISMIP7_ESM"] == "CESM2-WACCM"
    # no reference: the historical carries none and one cannot be built
    # from its evolved endpoint (simulation.setup_model refuses)
    assert c["ISMIP7_APPARENT_MB"] == "0" and c["ISMIP7_RUN_TAG"] == "out"
    assert c["ISMIP7_FRACTURE"] == "none" and "ISMIP7_T_START" not in c
    assert "ISMIP7_GEOMETRY_BACKDATE" not in c
    p = by_name(subs, "rt_out_ssp585_cesm_waccm")["env"]
    assert p["ISMIP7_EXPERIMENT"] == "ssp585_cesm_waccm"
    assert p["ISMIP7_APPARENT_MB"] == "0" and p["ISMIP7_FRACTURE"] == "mask_front"

    finish_forward(results / "hist_mri_esm2_0_out_32000_final.h5", 2015)
    routine(sandbox, "next")
    names = sorted(s["argv"].split(" -J ")[1].split()[0] for s in submits(sandbox))
    assert names == ["rt_out_control_mri_esm2", "rt_out_ssp126_mri_esm2",
                     "rt_out_ssp370_mri_esm2", "rt_out_ssp585_mri_esm2"]

    rc, log = routine(sandbox, "next")
    assert rc == 0 and submits(sandbox) == []
    rc, log = routine(sandbox, "status")
    assert rc == 0, log
    assert "r_friction2          done (20 its" in log
    assert "hist_cesm_waccm      done (t=2015.0)" in log
    assert "ssp585_mri_esm2      started, no checkpoint yet" in log


def test_init_refuses_to_overwrite_a_recorded_routine(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, _ = routine(sandbox, "init")
    assert rc == 0
    submits(sandbox)
    rc, log = routine(sandbox, "init")
    assert rc != 0 and "exists" in log
    assert submits(sandbox) == []


def test_slurm_mode_needs_its_resources(sandbox):
    rc, log = routine(sandbox, "init", env={"ISMIP7_ROUTINE_SBATCH_FWD": ""})
    assert rc != 0 and "ISMIP7_ROUTINE_SBATCH_FWD" in log


def test_forwards_alone_run_from_the_map_given(sandbox):
    r"""Without an inversion stage the forwards start from
    ISMIP7_ROUTINE_WARM_START: the routine as a plain forward suite."""
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={
        "ISMIP7_ROUTINE_STAGES": "hist control", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM"})
    assert rc == 0, log
    assert os.readlink(sandbox / "out" / "map_forward.h5") == str(sandbox / "warm.h5")
    (h,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in h["argv"]
    assert h["env"]["ISMIP7_INVERSION"] == str(sandbox / "out" / "map_forward.h5")


def test_resume_submits_a_started_stage_again(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    routine(sandbox, "init")
    submits(sandbox)
    rc, log = routine(sandbox, "resume", "shelf1")
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_shelf1" in again["argv"]
    assert again["env"]["ISMIP7_CHAIN_ITER_BUDGET"] == "40"
    rc, log = routine(sandbox, "resume", "friction1")
    assert rc != 0 and "not ready" in log


def test_a_2003_historical_keeps_the_default_backdate(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={
        "ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
        "ISMIP7_ROUTINE_HIST_START": "2003"})
    assert rc == 0, log
    (h,) = submits(sandbox)
    assert h["env"]["ISMIP7_T_START"] == "2003"
    assert "ISMIP7_GEOMETRY_BACKDATE" not in h["env"]


@pytest.mark.parametrize("precond", ["none", "mass", ""])
def test_init_refuses_inversion_stages_off_the_tao_path(sandbox, precond):
    (sandbox / "warm.h5").write_text("warm\n")
    env = {"ISMIP7_GRAD_PRECOND": precond}
    rc, log = routine(sandbox, "init", env=env)
    assert rc != 0 and "ISMIP7_GRAD_PRECOND" in log, log
    assert submits(sandbox) == []
    assert not (sandbox / "out" / "routine.env").exists(), "a refused init records nothing"


def test_forwards_alone_need_no_tao_metric(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={
        "ISMIP7_GRAD_PRECOND": "none", "ISMIP7_ROUTINE_STAGES": "hist",
        "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM"})
    assert rc == 0, log
    assert len(submits(sandbox)) == 1


def test_a_lock_whose_holder_is_gone_is_taken(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    routine(sandbox, "init")
    submits(sandbox)
    dead = subprocess.Popen(["true"])
    dead.wait()
    lock = sandbox / "out" / ".routine.lock"
    lock.mkdir()
    (lock / "holder").write_text(f"{dead.pid} {os.uname().nodename}\n")
    finish_inversion(sandbox / "out" / "map_shelf1.h5")
    rc, log = routine(sandbox, "next")
    assert rc == 0, log
    (s,) = submits(sandbox)
    assert "-J rt_out_friction1" in s["argv"]
    assert not lock.exists()


def test_a_live_lock_is_waited_for(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    routine(sandbox, "init")
    submits(sandbox)
    lock = sandbox / "out" / ".routine.lock"
    lock.mkdir()
    (lock / "holder").write_text(f"{os.getpid()} {os.uname().nodename}\n")
    finish_inversion(sandbox / "out" / "map_shelf1.h5")
    with pytest.raises(subprocess.TimeoutExpired):
        merged = base_env(sandbox)
        subprocess.run(
            ["bash", str(sandbox / "antarctica" / "scripts" / "initialization_routine.sh"),
             str(sandbox / "out"), "next"],
            env=merged, cwd=str(sandbox), capture_output=True, text=True, timeout=8)
    assert submits(sandbox) == []
    assert (lock / "holder").read_text().split()[0] == str(os.getpid())


def test_watch_resubmits_a_dead_chain_and_cancels_a_silent_link(sandbox):
    r"""A started stage with no job left in the queue is resubmitted (a link
    killed by the scheduler queues no successor); a running forward link whose
    log has been silent longer than the stall limit is cancelled and its node
    excluded from the stage's next submission."""
    import time
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    queue = sandbox / "queue.txt"
    (bin_dir / "squeue").write_text(
        "#!/bin/bash\n"
        'name=""; while [ $# -gt 0 ]; do [ "$1" = -n ] && name=$2; shift; done\n'
        '[ -f "$QUEUE_FILE" ] && grep "^$name " "$QUEUE_FILE" | cut -d" " -f2-\n'
        "exit 0\n")
    (bin_dir / "scancel").write_text('#!/bin/bash\necho "scancel $*" >> "$SBATCH_CALLS"\n')
    for name in ("squeue", "scancel"):
        (bin_dir / name).chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist control", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "QUEUE_FILE": str(queue), "ISMIP7_ROUTINE_WATCH_ONCE": "1",
           "ISMIP7_ROUTINE_STALL_MIN": "15"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in h["argv"]

    # the chain died: nothing queued under the stage's name
    queue.write_text("")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert "--exclude" not in again["argv"]
    assert "no job in the queue; resubmitting" in log

    # a running link whose log went silent 20 minutes ago
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    silent = logs / "ismip7_fwd_4242.out"
    silent.write_text("t=2003.0\n")
    old = time.time() - 20 * 60
    os.utime(silent, (old, old))
    queue.write_text("rt_out_hist_cesm_waccm 4242 RUNNING bb5u20c1 0:30:00\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert "silent for 20 min" in log
    calls = (sandbox / "sbatch_calls.txt").read_text()
    assert "scancel 4242" in calls
    assert (sandbox / "out" / "exclude_nodes").read_text().strip() == "bb5u20c1"
    # the cancelled job is gone on the next pass: the stage is resubmitted
    # without that node
    queue.write_text("")
    (sandbox / "sbatch_calls.txt").unlink()
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "--exclude=bb5u20c1" in again["argv"]

    # a running link with a fresh log is left alone
    silent.write_text("t=2004.0\n")
    queue.write_text("rt_out_hist_cesm_waccm 4243 RUNNING bb5u21c1 0:10:00\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == []
    assert "scancel 4243" not in (sandbox / "sbatch_calls.txt").read_text() if (sandbox / "sbatch_calls.txt").exists() else True


def test_watch_does_nothing_when_the_scheduler_cannot_be_reached(sandbox):
    r"""An unreachable controller answers nothing; that is not an empty queue.
    The 14 October 2026 outage resubmitted every stage on top of its live chain."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    (bin_dir / "squeue").write_text('#!/bin/bash\necho "squeue: error: Unable to contact slurm controller" >&2\nexit 1\n')
    (bin_dir / "squeue").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert "squeue failed; skipping" in log
    assert submits(sandbox) == [], "nothing may be resubmitted on a scheduler failure"



def test_watch_parks_a_stage_whose_driver_failed(sandbox):
    r"""A link that exits with a driver error ends its chain like a scheduler
    kill does, but a resubmission only repeats the error: on 15 October 2026
    each control and projection of two routines was resubmitted about 60
    times, two minutes apart, on an a_ref refusal. The watchdog parks such a
    stage with the error, status shows it, and resume clears it."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    queue = sandbox / "queue.txt"
    (bin_dir / "squeue").write_text(
        "#!/bin/bash\n"
        'name=""; while [ $# -gt 0 ]; do [ "$1" = -n ] && name=$2; shift; done\n'
        '[ -f "$QUEUE_FILE" ] && grep "^$name " "$QUEUE_FILE" | cut -d" " -f2-\n'
        "exit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "QUEUE_FILE": str(queue), "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    jid = (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    (logs / f"ismip7_fwd_{jid}.out").write_text(
        "=== forward exit 1 ===\nchain: stopping, driver exited 1\n")
    (logs / f"ismip7_fwd_{jid}.err").write_text(
        "Traceback\nRuntimeError: ISMIP7_APPARENT_MB is set but restart checkpoint x has no a_ref_mb\n")
    queue.write_text("")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == [], "a driver failure is not resubmitted"
    assert "stopped on a driver failure" in log
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert "no a_ref_mb" in blocked
    rc, log = routine(sandbox, "watch", env=env)
    assert submits(sandbox) == [], "a parked stage stays parked"
    rc, log = routine(sandbox, "status", env=env)
    assert "BLOCKED" in log and "no a_ref_mb" in log
    rc, log = routine(sandbox, "resume", "hist_cesm_waccm", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists()


def _stalled_link(sandbox, resumed, stalled):
    r"""A stage whose newest forward link resumed at ``resumed`` and stalled
    at ``stalled``, as projection.sbatch and the driver log it."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    (bin_dir / "squeue").write_text("#!/bin/bash\nexit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    jid = (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    (logs / f"ismip7_fwd_{jid}.out").write_text(
        f"  Restart: evolved geometry + frozen anchors loaded (t_yr={resumed}, friction=budd)\n"
        "  Step 1: rescue ladder + subcycles exhausted, saving and stopping\n"
        "=== 2026-10-16T12:11:54-05:00  forward exit 0 ===\n"
        f"chain: the solver stalled (rescue ladder + subcycles exhausted) at {stalled}.\n"
        "chain: NOT resubmitting. A stall is not a run that needs more time;\n")
    return env, jid


def test_watch_parks_a_stage_that_stalled_where_it_resumed(sandbox):
    r"""A forward whose solver stalled saves its state, exits 0 and does not
    chain. The watchdog used to read that as a link the scheduler took and
    resubmit it; a link that stalls where it resumed makes no progress, and on
    16 October 2026 the seq0 ssp126 projection was relaunched link after link,
    each stalling in its first step at 2019.95. Parked, with the year."""
    env, jid = _stalled_link(sandbox, "2019.95", "2019.95")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    assert submits(sandbox) == [], "a stall without progress is not resubmitted"
    blocked = (sandbox / "out" / "hist_cesm_waccm.blocked").read_text()
    assert blocked.strip() == f"job {jid}: the solver stalled at 2019.95, where the link resumed"
    rc, log = routine(sandbox, "resume", "hist_cesm_waccm", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]


def test_watch_resubmits_a_stage_that_advanced_before_it_stalled(sandbox):
    r"""A link that advanced and then stalled is resubmitted: the next link
    starts the substep controller afresh at the stall, and the seq0 ssp585
    projection's restart passed the step its link (2015.0 to 2021.25) had
    stalled in."""
    env, _ = _stalled_link(sandbox, "2015.0", "2021.25")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]
    assert not (sandbox / "out" / "hist_cesm_waccm.blocked").exists()


def test_a_scheduler_kill_is_still_resubmitted(sandbox):
    r"""A link the scheduler took (no driver failure in its log) is resubmitted."""
    (sandbox / "warm.h5").write_text("warm\n")
    bin_dir = sandbox / "bin"
    (bin_dir / "squeue").write_text("#!/bin/bash\nexit 0\n")
    (bin_dir / "squeue").chmod(0o755)
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_WATCH_ONCE": "1"}
    routine(sandbox, "init", env=env)
    submits(sandbox)
    jid = (sandbox / "out" / "hist_cesm_waccm.jobs").read_text().split()[0]
    logs = sandbox / "logs"
    logs.mkdir(exist_ok=True)
    (logs / f"ismip7_fwd_{jid}.out").write_text("  t=2005.3 ...\n")
    rc, log = routine(sandbox, "watch", env=env)
    assert rc == 0, log
    (again,) = submits(sandbox)
    assert "-J rt_out_hist_cesm_waccm" in again["argv"]


def test_the_forwards_run_the_routine_s_calving_law(sandbox):
    r"""The forwards take the routine's calving law, its parameters and its
    advance rule, recorded at init; the defaults are the submission front
    (hfb at 150 kPa, retreat-only), the collapse mask stays with the
    projections, and the legacy pinned front stays off under a law."""
    (sandbox / "warm.h5").write_text("warm\n")
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_CALVING": "vonmises",
           "ISMIP7_ROUTINE_CALVING_PARAMS": "sigma_max_fl=0.275,sigma_max_gr=1.0",
           "ISMIP7_ROUTINE_FRONT_ADVANCE": "free"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    e = h["env"]
    assert e["ISMIP7_CALVING"] == "vonmises"
    assert e["ISMIP7_CALVING_PARAMS"] == "sigma_max_fl=0.275,sigma_max_gr=1.0"
    assert e["ISMIP7_FRONT_ADVANCE"] == "free"
    assert e["ISMIP7_FRACTURE"] == "none"
    assert "ISMIP7_FIXED_FRONT" not in e
    env_text = (sandbox / "out" / "routine.env").read_text()
    assert "ISMIP7_ROUTINE_CALVING=vonmises" in env_text


def test_a_law_without_parameters_runs_on_its_defaults(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM",
           "ISMIP7_ROUTINE_CALVING": "nocalving", "ISMIP7_ROUTINE_CALVING_PARAMS": "",
           "ISMIP7_ROUTINE_FRONT_ADVANCE": "free"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    assert h["env"]["ISMIP7_CALVING"] == "nocalving"
    assert "ISMIP7_CALVING_PARAMS" not in h["env"]


def test_a_routine_recorded_before_the_calving_knobs_runs_the_submission_front(sandbox):
    r"""The routines started on 13-14 October recorded no calving knobs; their
    forwards ran the submission front (ismip7_front_routine.sh) and must keep
    it, not be handed an empty law."""
    (sandbox / "warm.h5").write_text("warm\n")
    env = {"ISMIP7_ROUTINE_STAGES": "hist", "ISMIP7_ROUTINE_ESMS": "CESM2-WACCM"}
    rc, log = routine(sandbox, "init", env=env)
    assert rc == 0, log
    submits(sandbox)
    envfile = sandbox / "out" / "routine.env"
    envfile.write_text("".join(
        line for line in envfile.read_text().splitlines(keepends=True)
        if not line.startswith(("export ISMIP7_ROUTINE_CALVING",
                                "export ISMIP7_ROUTINE_FRONT_ADVANCE"))))
    rc, log = routine(sandbox, "resume", "hist_cesm_waccm", env=env)
    assert rc == 0, log
    (h,) = submits(sandbox)
    e = h["env"]
    assert e["ISMIP7_CALVING"] == "hfb"
    assert e["ISMIP7_CALVING_PARAMS"] == "sigma_max=0.15"
    assert e["ISMIP7_FRONT_ADVANCE"] == "none"
    assert "ISMIP7_FIXED_FRONT" not in e


def test_an_unknown_advance_rule_is_refused_at_init(sandbox):
    (sandbox / "warm.h5").write_text("warm\n")
    rc, log = routine(sandbox, "init", env={"ISMIP7_ROUTINE_FRONT_ADVANCE": "sometimes"})
    assert rc != 0 and "ISMIP7_ROUTINE_FRONT_ADVANCE" in log
    assert submits(sandbox) == []
