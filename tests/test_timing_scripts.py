r"""The timing campaign's job scripts, run against a Slurm shim.

They used to carry one cluster's account, modules and venv. They now source
`site_core.sh` and nothing else, so what is pinned here is what that must and
must not do to a lane:

- the site's Firedrake is activated and ranks start through `ismip7_mpirun`;
- no model default leaks in. `site_env.sh` exports ISMIP7_MESH, which
  `run_timing.py` refuses on a matrix lane, and defaults ISMIP7_FRICTION to a
  law the campaign does not run. A lane's model settings are its exports alone;
- the kernel cache persists between jobs, because a lane's seconds_per_step has
  no warm-up excluded;
- the status file still walks running -> finished/failed as the manager expects.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BR = "scripts/batch_runners"

DRIVER = r'''
import os, sys
with open(os.environ["FAKE_ENV_DUMP"], "a") as fh:
    fh.write("ARGV " + " ".join(sys.argv[1:]) + "\n")
    for key in sorted(os.environ):
        if key.startswith(("ISMIP7_", "PYOP2_", "FIREDRAKE_", "OMP_", "OPENBLAS_", "FAKE_VENV")):
            fh.write(f"{key}={os.environ[key]}\n")
status = os.environ.get("FAKE_WRITE_STATUS")
if status:
    open(os.environ["ISMIP7_TIMING_STATUS"], "w").write(status + "\n")
sys.exit(int(os.environ.get("FAKE_RC", "0")))
'''

STUBS = {
    # `srun -n N python -u script args...`: record N, hand the rest to the driver.
    "srun": ('#!/bin/bash\n'
             'echo "SRUN $*" >> "$FAKE_ENV_DUMP"\n'
             'shift 2; shift 2\n'
             'exec "$FAKE_PYTHON" "$FAKE_DRIVER" "$@"\n'),
    "scontrol": "#!/bin/bash\nexit 0\n",
    "module": "#!/bin/bash\nexit 0\n",
}


@pytest.fixture
def sandbox(tmp_path):
    r"""A submit directory that looks like antarctica/ to the job script."""
    root = tmp_path / "repo" / "antarctica"
    root.mkdir(parents=True)
    (root / "scripts").symlink_to(REPO / "antarctica" / "scripts")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in STUBS.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    (tmp_path / "driver.py").write_text(DRIVER)
    (tmp_path / "activate").write_text("export FAKE_VENV_ACTIVE=1\n")
    return tmp_path


def run_script(sandbox, name, **env):
    root = sandbox / "repo" / "antarctica"
    base = {
        "PATH": f"{sandbox / 'bin'}:{os.environ['PATH']}",
        "HOME": str(sandbox),
        "SLURM_JOB_ID": "777",
        "SLURM_NTASKS": "16",
        "SLURM_SUBMIT_DIR": str(root),
        "ISMIP7_SITE": "local",
        "ISMIP7_LOCAL_ENV": os.devnull,
        "ISMIP7_FIREDRAKE": str(sandbox / "activate"),
        "ISMIP7_REPO": str(sandbox / "repo"),
        "FAKE_PYTHON": sys.executable,
        "FAKE_DRIVER": str(sandbox / "driver.py"),
        "FAKE_ENV_DUMP": str(sandbox / "env_dump.txt"),
    }
    base.update(env)
    proc = subprocess.run(["bash", str(root / BR / name)], env=base, cwd=str(sandbox),
                          capture_output=True, text=True)
    dump = sandbox / "env_dump.txt"
    return proc, (dump.read_text().splitlines() if dump.exists() else [])


def test_a_lane_runs_in_the_site_s_environment_with_only_its_own_model_settings(sandbox):
    status = sandbox / "status.txt"
    proc, seen = run_script(sandbox, "timing_transient.script",
                            ISMIP7_TIMING_STATUS=str(status), ISMIP7_LC="500")
    assert proc.returncode == 0, proc.stderr
    assert "SRUN -n 16 python -u scripts/run_timing.py" in seen
    assert "FAKE_VENV_ACTIVE=1" in seen
    assert "OMP_NUM_THREADS=1" in seen and "OPENBLAS_NUM_THREADS=1" in seen
    assert "ISMIP7_LC=500" in seen
    leaked = [line for line in seen if line.startswith((
        "ISMIP7_MESH=", "ISMIP7_FRICTION=", "ISMIP7_LC_COARSE=", "ISMIP7_MAP_DEFAULT=",
        "ISMIP7_GEOMETRY_SPACE=", "ISMIP7_N_FLOW=", "ISMIP7_DATA_ROOT="))]
    assert leaked == []
    assert "site    local" in proc.stdout and "friction=" not in proc.stdout
    assert status.read_text().startswith("finished exit_code=0 job_id=777 ")


def test_the_kernel_cache_persists_between_lanes(sandbox):
    r"""By default Firedrake's own location (no PYOP2_CACHE_DIR, and no empty
    per-job directory left behind); a named one when the site gives it."""
    proc, seen = run_script(sandbox, "timing_transient.script")
    assert proc.returncode == 0, proc.stderr
    assert not any(line.startswith("PYOP2_CACHE_DIR=") for line in seen)
    assert not (sandbox / ".pyop2_cache" / "777").exists()
    jit = sandbox / "jit"
    proc, seen = run_script(sandbox, "timing_transient.script",
                            ISMIP7_TIMING_JIT_CACHE=str(jit))
    assert proc.returncode == 0, proc.stderr
    assert f"PYOP2_CACHE_DIR={jit}/pyop2" in seen
    assert f"FIREDRAKE_TSFC_KERNEL_CACHE_DIR={jit}/tsfc" in seen
    assert (jit / "pyop2").is_dir()


def test_a_killed_lane_is_stamped_failed_and_a_specific_stamp_is_kept(sandbox):
    status = sandbox / "status.txt"
    proc, _ = run_script(sandbox, "timing_transient.script",
                         ISMIP7_TIMING_STATUS=str(status), FAKE_RC="137")
    assert proc.returncode == 137
    assert status.read_text().startswith(
        "failed category=external_termination exit_code=137 job_id=777 ")
    proc, _ = run_script(sandbox, "timing_transient.script",
                         ISMIP7_TIMING_STATUS=str(status), FAKE_RC="1",
                         FAKE_WRITE_STATUS="failed category=tripwire phase=step_3")
    assert proc.returncode == 1
    assert status.read_text() == "failed category=tripwire phase=step_3\n"


def test_prepare_repacks_on_one_rank_through_the_same_launcher(sandbox):
    cache = sandbox / "cache"
    cache.mkdir()
    status = cache / "status.txt"
    final, manifest = cache / "state.h5", cache / "state.json"
    final.write_text("x"); manifest.write_text("{}")
    proc, seen = run_script(
        sandbox, "timing_prepare.script",
        ISMIP7_TIMING_CACHE_RAW=str(cache / "raw.h5"), ISMIP7_TIMING_CACHE=str(final),
        ISMIP7_TIMING_CACHE_MANIFEST=str(manifest), ISMIP7_TIMING_CACHE_STATUS=str(status),
        ISMIP7_TIMING_CACHE_PRISTINE=str(cache / "state.prepare.h5"))
    assert proc.returncode == 0, proc.stderr
    launches = [line for line in seen if line.startswith("SRUN ")]
    assert launches[0] == "SRUN -n 16 python -u scripts/prepare_timing_cache.py"
    assert launches[1].startswith("SRUN -n 1 python -u scripts/redistribute_checkpoint.py --input ")
    assert status.read_text().startswith("finished phase=published exit_code=0 ")
    assert (cache / "state.prepare.h5").exists()


@pytest.mark.parametrize("name", sorted(
    p.name for p in (REPO / "antarctica" / BR).glob("*.script")))
def test_no_job_script_names_a_cluster_a_person_or_a_resource(name):
    text = (REPO / "antarctica" / BR / name).read_text()
    for needle in ("/N/u/", "r00905", "@iu.edu", "module load", "srun "):
        assert needle not in text, f"{name} still carries {needle!r}"
    directives = [line for line in text.splitlines() if line.startswith("#SBATCH")]
    assert len(directives) == 2
    assert all(line.split()[1] in ("-o", "-e") for line in directives), directives
    assert ". scripts/batch_runners/site_core.sh" in text
    assert ". scripts/batch_runners/site_env.sh" not in text
