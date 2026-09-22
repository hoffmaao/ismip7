r"""A site whose Firedrake is a container image rather than a venv.

UChicago Midway runs from an Apptainer image. Setting ISMIP7_CONTAINER has to be
the whole difference: the same job scripts, with `python` resolving to a shim
that runs in the image, and ranks started by the image's own mpiexec inside one
`exec` (every job is one node, so the MPI in the image never has to agree with
the host's Slurm about PMI). These run the shipped timing script and the shim
against a fake `apptainer` that records how it was called and then runs the
command on the host.
"""
import os
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BR = REPO / "antarctica" / "scripts" / "batch_runners"

# apptainer exec [--bind X]... IMAGE cmd...: log it, drop through to cmd.
FAKE_APPTAINER = r'''#!/bin/bash
echo "APPTAINER $*" >> "$FAKE_LOG"
[ "$1" = exec ] || exit 64
shift
while [ "${1#--}" != "$1" ]; do
    case "$1" in --bind) shift 2 ;; *) shift ;; esac
done
image="$1"; shift
[ -r "$image" ] || exit 65
exec "$@"
'''
FAKE_MPIEXEC = r'''#!/bin/bash
echo "MPIEXEC $*" >> "$FAKE_LOG"
while [ "$1" != python3 ]; do shift; done
exec "$@"
'''


@pytest.fixture
def sandbox(tmp_path):
    root = tmp_path / "repo" / "antarctica"
    root.mkdir(parents=True)
    (root / "scripts").symlink_to(REPO / "antarctica" / "scripts")
    (tmp_path / "work").mkdir()
    (tmp_path / "image.sif").write_text("not really an image\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stubs = {
        "apptainer": FAKE_APPTAINER,
        "mpiexec": FAKE_MPIEXEC,
        # The image's python3: log the argv and the environment that reached it.
        "python3": ('#!/bin/bash\necho "PYTHON3 $*" >> "$FAKE_LOG"\n'
                    'env | grep -E "^(ISMIP7_LC|OMP_NUM_THREADS|PYOP2_CACHE_DIR)=" >> "$FAKE_LOG"\n'
                    'exit "${FAKE_RC:-0}"\n'),
        "srun": '#!/bin/bash\necho "SRUN $*" >> "$FAKE_LOG"\nexit 99\n',
        "scontrol": "#!/bin/bash\nexit 0\n",
        "module": "#!/bin/bash\nexit 0\n",
    }
    for name, body in stubs.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    return tmp_path


def env_for(sandbox, **extra):
    base = {
        "PATH": f"{sandbox / 'bin'}:/usr/bin:/bin",
        "HOME": str(sandbox),
        "SLURM_JOB_ID": "888",
        "SLURM_NTASKS": "16",
        "SLURM_SUBMIT_DIR": str(sandbox / "repo" / "antarctica"),
        "ISMIP7_LOCAL_ENV": os.devnull,
        "ISMIP7_SITE": "uchicago_midway",
        "ISMIP7_PART_LONG": "caslake", "ISMIP7_PART_SHORT": "caslake",
        "ISMIP7_PART_DEBUG": "caslake",
        "ISMIP7_REPO": str(sandbox / "repo"),
        "ISMIP7_WORK": str(sandbox / "work"),
        "ISMIP7_CONTAINER": str(sandbox / "image.sif"),
        "ISMIP7_TIMING_JIT_CACHE": str(sandbox / "work" / "jit"),
        "FAKE_LOG": str(sandbox / "log.txt"),
    }
    base.update(extra)
    return base


def log(sandbox):
    path = sandbox / "log.txt"
    return path.read_text().splitlines() if path.exists() else []


def test_a_timing_lane_runs_in_the_image_with_no_venv_and_no_srun(sandbox):
    status = sandbox / "status.txt"
    scratch = sandbox / "scratch"
    proc = subprocess.run(
        ["bash", str(sandbox / "repo/antarctica/scripts/batch_runners/timing_transient.script")],
        env=env_for(sandbox, ISMIP7_TIMING_STATUS=str(status), ISMIP7_LC="500",
                    SCRATCH=str(scratch)),
        cwd=str(sandbox), capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    seen = log(sandbox)
    assert not any(line.startswith("SRUN") for line in seen)
    launch = next(line for line in seen if line.startswith("APPTAINER"))
    image = sandbox / "image.sif"
    assert launch.endswith(f"{image} mpiexec -n 16 python3 -u scripts/run_timing.py")
    # The per-job loopy cache is on SCRATCH, outside the $HOME the runtime binds
    # itself, so it needs a bind of its own exactly as PyOP2's shared one does.
    for bound in (sandbox / "repo", sandbox / "work", sandbox / "work" / "jit",
                  scratch / ".pyop2_cache" / "xdg" / "888"):
        assert f"--bind {bound}" in launch, launch
    assert "--cleanenv" not in launch
    assert "MPIEXEC -n 16 python3 -u scripts/run_timing.py" in seen
    assert "ISMIP7_LC=500" in seen and "OMP_NUM_THREADS=1" in seen
    assert f"PYOP2_CACHE_DIR={sandbox}/work/jit/pyop2" in seen
    assert status.read_text().startswith("finished exit_code=0 ")


def test_a_failure_in_the_image_is_the_job_s_failure(sandbox):
    status = sandbox / "status.txt"
    proc = subprocess.run(
        ["bash", str(sandbox / "repo/antarctica/scripts/batch_runners/timing_transient.script")],
        env=env_for(sandbox, ISMIP7_TIMING_STATUS=str(status), FAKE_RC="9"),
        cwd=str(sandbox), capture_output=True, text=True)
    assert proc.returncode == 9
    assert status.read_text().startswith("failed category=external_termination exit_code=9 ")


def test_serial_python_lines_reach_the_image_through_the_shim(sandbox):
    script = ('. "%s/site_core.sh"; ismip7_activate; ismip7_activate; '
              'echo "PATH=$PATH"; python -c "print(1)" an-arg' % BR)
    proc = subprocess.run(["bash", "-c", script], env=env_for(sandbox),
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    seen = log(sandbox)
    assert seen[0].startswith("APPTAINER exec --bind ")
    assert seen[0].endswith(f'{sandbox / "image.sif"} python3 -c print(1) an-arg')
    # Activating twice (a chain successor inherits PATH) adds the shim once.
    path = proc.stdout.split("PATH=", 1)[1].splitlines()[0]
    assert path.split(":").count(str(BR / "container_bin")) == 1
    assert path.split(":")[0] == str(BR / "container_bin")


def test_the_launcher_inside_the_image_can_be_given_flags(sandbox):
    script = '. "%s/site_core.sh"; ismip7_activate; ismip7_mpirun 1 python -u x.py --input a' % BR
    proc = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True,
        env=env_for(sandbox, ISMIP7_CONTAINER_MPIEXEC="mpiexec --mca plm isolated",
                    ISMIP7_CONTAINER_ARGS="--bind /project2"))
    assert proc.returncode == 0, proc.stderr
    launch = log(sandbox)[0]
    assert "--bind /project2" in launch
    assert launch.endswith("mpiexec --mca plm isolated -n 1 python3 -u x.py --input a")


def test_a_missing_image_or_runtime_is_said_plainly(sandbox):
    script = '. "%s/site_core.sh"; ismip7_activate' % BR
    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                          env=env_for(sandbox, ISMIP7_CONTAINER=str(sandbox / "nope.sif")))
    assert proc.returncode == 2 and "ISMIP7_CONTAINER is unreadable" in proc.stderr
    (sandbox / "bin" / "apptainer").unlink()
    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                          env=env_for(sandbox))
    assert proc.returncode == 2 and "'apptainer' is not on PATH" in proc.stderr


def test_submit_asks_a_container_site_for_no_venv(sandbox):
    proc = subprocess.run(
        ["bash", str(BR / "submit.sh"), "script", "scripts/batch_runners/timing_transient.script",
         "--cd", "antarctica", "--dry-run"],
        env=env_for(sandbox), capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert "site uchicago_midway: env ISMIP7_SITE=uchicago_midway " in proc.stderr
    assert " -p caslake " in proc.stderr
    # Without an image the site is still the stub it ships as.
    env = env_for(sandbox)
    for name in ("ISMIP7_CONTAINER", "ISMIP7_PART_LONG", "ISMIP7_PART_SHORT", "ISMIP7_PART_DEBUG"):
        del env[name]
    proc = subprocess.run(["bash", str(BR / "submit.sh"), "smoke"], env=env,
                          capture_output=True, text=True)
    assert proc.returncode == 2
    assert "ISMIP7_FIREDRAKE" in proc.stderr and "ISMIP7_PART_DEBUG" in proc.stderr


def test_the_shim_refuses_to_run_outside_a_container_site(sandbox):
    env = env_for(sandbox)
    del env["ISMIP7_CONTAINER"]
    proc = subprocess.run([str(BR / "container_bin" / "python"), "-c", "1"], env=env,
                          capture_output=True, text=True)
    assert proc.returncode == 2 and "ISMIP7_CONTAINER is not set" in proc.stderr


def test_an_image_that_inherits_the_host_path_cannot_recurse(sandbox):
    r"""Seen while writing these tests: a runtime that hands the host PATH to
    the image finds the shim again as `python` and starts a container per call
    without end. The shim marks the environment and refuses the second entry."""
    proc = subprocess.run([str(BR / "container_bin" / "python"), "-c", "1"],
                          env=env_for(sandbox, ISMIP7_IN_CONTAINER="1"),
                          capture_output=True, text=True, timeout=30)
    assert proc.returncode == 2 and "reached from inside the image" in proc.stderr
    assert log(sandbox) == []
    assert not (BR / "container_bin" / "python3").exists()
