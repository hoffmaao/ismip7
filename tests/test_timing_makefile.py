r"""antarctica/Makefile's own submissions, against a fake sbatch.

`transient-direct` (the qualification, debug and reference lanes), `meshes` and
`redistribute` submit through `submit.sh script` like the campaign manager
does. The recipe is shell inside make, so this runs it for real in a sandbox
whose results/ is a temporary directory: the site's account and partition
class reach sbatch, the job starts from antarctica/, the status file is
stamped with the job id, and a lane this site's nodes cannot hold is recorded
not_runnable without failing the target.
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(shutil.which("make") is None, reason="needs make")

FAKE_SBATCH = '''#!/bin/bash
printf 'CWD: %s\\n' "$PWD" >> "$SBATCH_CALLS"
printf 'ARG: %s\\n' "$@" >> "$SBATCH_CALLS"
echo "5151;cluster"
'''


@pytest.fixture
def sandbox(tmp_path):
    root = tmp_path / "repo" / "antarctica"
    root.mkdir(parents=True)
    (root / "Makefile").symlink_to(REPO / "antarctica" / "Makefile")
    (root / "scripts").symlink_to(REPO / "antarctica" / "scripts")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "sbatch").write_text(FAKE_SBATCH)
    (bin_dir / "sbatch").chmod(0o755)
    (tmp_path / "input").write_text("stands in for a mesh, a sidecar and a MAP\n")
    return tmp_path


def make(sandbox, target, *variables, **env):
    root = sandbox / "repo" / "antarctica"
    base = {
        "PATH": f"{sandbox / 'bin'}:{os.environ['PATH']}",
        "HOME": str(sandbox),
        "ISMIP7_LOCAL_ENV": os.devnull,
        "ISMIP7_FIREDRAKE": str(sandbox / "input"),
        "SBATCH_CALLS": str(sandbox / "sbatch_calls.txt"),
    }
    base.update(env)
    proc = subprocess.run(["make", "-s", target, *variables], cwd=str(root), env=base,
                          capture_output=True, text=True)
    calls = sandbox / "sbatch_calls.txt"
    return proc, (calls.read_text().splitlines() if calls.exists() else [])


def lane(sandbox, cores="16"):
    stand_in = str(sandbox / "input")
    return ["LCS=2500", "RATIOS=10", f"CORES={cores}", "TIMING_KIND=debug",
            "TIMING_TAG=maketest", f"TIMING_MESH={stand_in}",
            f"TIMING_BNDIDS={stand_in}", f"TIMING_INVERSION={stand_in}", "FORCE_TIMING=1"]


def status(sandbox, cores="16"):
    path = sandbox / "repo/antarctica/results/timing" / f"status_maketest_2500_25000_{cores}.txt"
    return path.read_text()


def test_a_direct_lane_goes_through_the_site_file(sandbox):
    proc, seen = make(sandbox, "transient-direct", *lane(sandbox), ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    root = sandbox / "repo" / "antarctica"
    assert seen[0] == f"CWD: {root}"
    assert seen[-1] == "ARG: scripts/batch_runners/timing_transient.script"
    for arg in ("-A", "r00905", "general", "--ntasks-per-node=16", "--mem=64G",
                "--time=12:00:00", "timing_2500_25000_16"):
        assert f"ARG: {arg}" in seen, seen
    export = next(line for line in seen if line.startswith("ARG: --export="))
    assert f"ISMIP7_REPO={sandbox / 'repo'}" in export
    assert f"ISMIP7_TIMING_STATUS={root}/results/timing/status_maketest_2500_25000_16.txt" in export
    assert "ISMIP7_LC=2500" in export and "ISMIP7_FRICTION=budd" in export
    assert status(sandbox).startswith("submitted job_id=5151 ")
    assert "site iu_quartz: sbatch" in proc.stderr


def test_the_debug_lanes_ask_for_the_site_s_debug_partition(sandbox):
    proc, seen = make(sandbox, "transient-direct", *lane(sandbox), "SLURM_QUEUE=debug",
                      ISMIP7_SITE="rice_nots")
    assert proc.returncode == 0, proc.stderr
    assert "ARG: scavenge" in seen and "ARG: -A" not in seen
    dry = subprocess.run(["make", "-n", "debug"], cwd=str(sandbox / "repo/antarctica"),
                         capture_output=True, text=True)
    assert '--queue "debug"' in dry.stdout


def test_a_partition_by_name_still_wins(sandbox):
    proc, seen = make(sandbox, "transient-direct", *lane(sandbox), "SLURM_PARTITION=gpu",
                      ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert "ARG: gpu" in seen and "ARG: general" not in seen


def test_a_lane_too_big_for_the_site_is_not_runnable_and_not_a_failure(sandbox):
    proc, seen = make(sandbox, "transient-direct", *lane(sandbox, "64"),
                      ISMIP7_SITE="rice_nots", ISMIP7_CORES_PER_NODE="40")
    assert proc.returncode == 0, proc.stderr
    assert seen == []
    assert status(sandbox, "64").startswith("not_runnable reason=exceeds_site_cores ")


def test_a_refused_submission_fails_the_target(sandbox):
    (sandbox / "bin" / "sbatch").write_text("#!/bin/bash\necho refused >&2\nexit 1\n")
    proc, _ = make(sandbox, "transient-direct", *lane(sandbox), ISMIP7_SITE="iu_quartz")
    assert proc.returncode != 0
    assert status(sandbox).startswith("submission_failed ")


def test_meshes_and_redistribute_wait_on_a_one_rank_job(sandbox):
    proc, seen = make(sandbox, "meshes", ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert seen[-1] == "ARG: scripts/batch_runners/timing_meshes.script"
    for arg in ("--wait", "--ntasks-per-node=1", "--mem=32G", "ant_timing_mesh", "r00905"):
        assert f"ARG: {arg}" in seen, seen
    export = next(line for line in seen if line.startswith("ARG: --export="))
    assert "TIMING_LCS=500:1000:2000:2500:5000" in export and "TIMING_BUFFER=20000" in export


def source_maps(sandbox):
    inv = sandbox / "inv"
    inv.mkdir(exist_ok=True)
    return inv / "map.parallel.h5", inv / "map.h5"


def test_a_repacked_map_whose_parallel_original_is_gone_is_a_complete_source(sandbox):
    r"""The invert job repacks its 16-rank MAP to one rank and deletes the
    original, so in production only the repack exists. `inversion` demanded the
    original, and `make timing` stopped there before any stage."""
    raw, packed = source_maps(sandbox)
    packed.write_text("one-rank repack\n")
    proc, seen = make(sandbox, "redistribute", f"SOURCE_TIMING_INVERSION={raw}",
                      f"TIMING_INVERSION={packed}", ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert seen == []
    assert "its parallel original is gone" in proc.stdout and "up to date" in proc.stdout


def test_a_parallel_original_with_no_repack_is_still_redistributed(sandbox):
    raw, packed = source_maps(sandbox)
    raw.write_text("sixteen-rank original\n")
    proc, seen = make(sandbox, "redistribute", f"SOURCE_TIMING_INVERSION={raw}",
                      f"TIMING_INVERSION={packed}", ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert seen[-1] == "ARG: scripts/batch_runners/timing_redistribute.script"
    export = next(line for line in seen if line.startswith("ARG: --export="))
    assert f"TIMING_REDISTRIBUTE_INPUT={raw}" in export
    assert f"TIMING_REDISTRIBUTE_OUTPUT={packed}" in export


def test_no_map_at_all_names_both_files(sandbox):
    raw, packed = source_maps(sandbox)
    proc, seen = make(sandbox, "redistribute", f"SOURCE_TIMING_INVERSION={raw}",
                      f"TIMING_INVERSION={packed}", ISMIP7_SITE="iu_quartz")
    assert proc.returncode != 0 and seen == []
    assert str(raw) in proc.stderr and str(packed) in proc.stderr
