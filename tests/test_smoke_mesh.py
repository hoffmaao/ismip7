r"""`smoke.sbatch` is the first job anyone runs at a new site, and no mesh is
tracked, so it has to find or make its own.

Seen on IU Quartz (job 10503057): the smoke named `antarctica_320000_32000.msh`,
a file only Rice has, left by a generator that did not yet tag the outline
buffer, and died in `Mesh()` five seconds in. Pinned here:

- with no mesh at all it builds one under the name the generator writes today,
  before any rank starts, and the ranks load that file;
- a mesh already there, under either name, is used and nothing is built;
- a mesh named through ISMIP7_MESH is the caller's: never built, never renamed;
- a failed build ends the job with the build's exit status, not an inversion
  that cannot open its mesh.
"""
import os
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SMOKE = "antarctica/scripts/batch_runners/smoke.sbatch"

STUBS = {
    # The serial `python` of the job: only the mesh generator is expected.
    "python": ('#!/bin/bash\n'
               'echo "PYTHON $*" >> "$FAKE_LOG"\n'
               '[ "${FAKE_MESH_RC:-0}" = 0 ] || exit "$FAKE_MESH_RC"\n'
               'touch "$FAKE_BUILDS"\n'),
    "srun": ('#!/bin/bash\n'
             'echo "SRUN $* MESH=$ISMIP7_MESH" >> "$FAKE_LOG"\n'
             '[ -f "$ISMIP7_MESH" ] || { echo "cannot open $ISMIP7_MESH" >&2; exit 65; }\n'),
    "scontrol": "#!/bin/bash\nexit 0\n",
    "module": "#!/bin/bash\nexit 0\n",
}


@pytest.fixture
def sandbox(tmp_path):
    repo = tmp_path / "repo"
    (repo / "antarctica" / "mesh").mkdir(parents=True)
    (repo / "antarctica" / "scripts").symlink_to(REPO / "antarctica" / "scripts")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in STUBS.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    (tmp_path / "activate").write_text("export FAKE_VENV_ACTIVE=1\n")
    return tmp_path


def run_smoke(sandbox, **env):
    repo = sandbox / "repo"
    mesh_dir = repo / "antarctica" / "mesh"
    base = {
        "PATH": f"{sandbox / 'bin'}:{os.environ['PATH']}",
        "HOME": str(sandbox),
        "SLURM_JOB_ID": "778",
        "SLURM_NTASKS": "4",
        "SLURM_SUBMIT_DIR": str(repo),
        "ISMIP7_SITE": "local",
        "ISMIP7_LOCAL_ENV": os.devnull,
        "ISMIP7_FIREDRAKE": str(sandbox / "activate"),
        "ISMIP7_REPO": str(repo),
        "FAKE_LOG": str(sandbox / "log.txt"),
        "FAKE_BUILDS": str(mesh_dir / "antarctica_320000_32000_buffered0.msh"),
    }
    base.update(env)
    proc = subprocess.run(["bash", str(repo / SMOKE)], env=base, cwd=str(sandbox),
                          capture_output=True, text=True, timeout=60)
    log = sandbox / "log.txt"
    return proc, (log.read_text().splitlines() if log.exists() else [])


def test_a_site_with_no_mesh_builds_the_smoke_mesh_first(sandbox):
    proc, seen = run_smoke(sandbox)
    assert proc.returncode == 0, proc.stderr
    built = sandbox / "repo" / "antarctica" / "mesh" / "antarctica_320000_32000_buffered0.msh"
    assert seen == [
        "PYTHON -u antarctica/scripts/mesh_antarctica.py --lc 32000 --lc-coarse 320000 --buffer-m 0",
        f"SRUN -n 4 python -u antarctica/scripts/inversion_icepack2.py MESH={built}",
    ]
    assert "building it" in proc.stdout and "smoke exit 0" in proc.stdout


@pytest.mark.parametrize("name", ["antarctica_320000_32000_buffered0.msh",
                                  "antarctica_320000_32000.msh"])
def test_a_mesh_already_there_is_used_as_it_is(sandbox, name):
    mesh = sandbox / "repo" / "antarctica" / "mesh" / name
    mesh.touch()
    proc, seen = run_smoke(sandbox)
    assert proc.returncode == 0, proc.stderr
    assert seen == [f"SRUN -n 4 python -u antarctica/scripts/inversion_icepack2.py MESH={mesh}"]


def test_a_mesh_the_caller_names_is_never_built(sandbox):
    proc, seen = run_smoke(sandbox, ISMIP7_MESH=str(sandbox / "elsewhere.msh"))
    assert proc.returncode == 65
    assert [line for line in seen if line.startswith("PYTHON")] == []
    assert seen[0].endswith(f"MESH={sandbox / 'elsewhere.msh'}")


def test_a_failed_build_ends_the_job_with_its_status(sandbox):
    proc, seen = run_smoke(sandbox, FAKE_MESH_RC="7")
    assert proc.returncode == 7
    assert [line for line in seen if line.startswith("SRUN")] == []
    assert "smoke mesh build exit 7" in proc.stdout
