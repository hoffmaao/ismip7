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
env | grep -E '^(ISMIP7|TIMING)_' | sort | sed 's/^/ENV: /' >> "$SBATCH_CALLS"
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


def environment(seen):
    r"""What submit.sh set in sbatch's environment, which a bare --export=ALL
    carries into the job."""
    return dict(line[len("ENV: "):].split("=", 1) for line in seen if line.startswith("ENV: "))


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
    assert "ARG: --export=ALL" in seen
    env = environment(seen)
    assert env["ISMIP7_REPO"] == str(sandbox / "repo")
    assert env["ISMIP7_TIMING_STATUS"] == f"{root}/results/timing/status_maketest_2500_25000_16.txt"
    assert env["ISMIP7_LC"] == "2500" and env["ISMIP7_FRICTION"] == "budd"
    assert status(sandbox).startswith("submitted job_id=5151 ")
    assert "site iu_quartz: env ISMIP7_SITE=iu_quartz " in proc.stderr


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
    env = environment(seen)
    assert env["TIMING_LCS"] == "500:1000:2000:2500:5000" and env["TIMING_BUFFER"] == "20000"


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
    env = environment(seen)
    assert env["TIMING_REDISTRIBUTE_INPUT"] == str(raw)
    assert env["TIMING_REDISTRIBUTE_OUTPUT"] == str(packed)


def test_no_map_at_all_names_both_files(sandbox):
    raw, packed = source_maps(sandbox)
    proc, seen = make(sandbox, "redistribute", f"SOURCE_TIMING_INVERSION={raw}",
                      f"TIMING_INVERSION={packed}", ISMIP7_SITE="iu_quartz")
    assert proc.returncode != 0 and seen == []
    assert str(raw) in proc.stderr and str(packed) in proc.stderr


@pytest.mark.parametrize("lc, dt, t_end", [("2000", "0.1", "2016"), ("5000", "0.125", "2016.25")])
def test_a_matrix_lane_s_step_scales_with_the_mesh_only_up_to_the_cap(sandbox, lc, dt, t_end):
    variables = [v for v in lane(sandbox) if not v.startswith(("LCS=", "TIMING_KIND="))]
    proc, seen = make(sandbox, "transient-direct", f"LCS={lc}", "TIMING_KIND=matrix",
                      *variables, ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    env = environment(seen)
    assert env["ISMIP7_DT"] == dt, seen
    assert env["ISMIP7_T_END"] == t_end, seen


def dry(sandbox, target, *variables):
    return subprocess.run(["make", "-n", target, *variables], capture_output=True, text=True,
                          cwd=str(sandbox / "repo" / "antarctica")).stdout


def test_the_campaign_solver_reaches_the_manager_and_the_qualification_gate(sandbox):
    r"""One variable names the solver for the whole path to a scout: the gate
    that qualifies it, the tag its records carry, and the manager's lanes."""
    assert '--solver "scpc_mumps"' in dry(sandbox, "timing-scout")
    for target in ("timing-scout", "timing-scale", "timing-probe"):
        assert '--solver "scpc_gamg"' in dry(sandbox, target, "TIMING_SOLVER=scpc_gamg"), target
    gate = dry(sandbox, "qualify-2step", "TIMING_SOLVER=scpc_gamg")
    assert "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_gamg" in gate
    assert "timing_qualification_2step_scpc_gamg_dg0_logvelnet_2500_64000_16.json" in gate
    assert "scpc_mumps" not in gate
    # Any mode solverconfig knows can still be put through the gate by name.
    assert "timing_qualification_2step_schur_gamg_" in dry(
        sandbox, "qualify-2step", "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=schur_gamg")


def test_make_matrix_renders_the_campaign_the_command_line_names(sandbox):
    r"""`make matrix` falls back to the campaign `make timing` last launched.
    A GAMG campaign started from `make timing-scout` never writes that file,
    so naming the solver has to win over it or the MUMPS matrix is rendered."""
    import sys
    (sandbox / "bin" / "python").write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    (sandbox / "bin" / "python").chmod(0o755)
    timing = sandbox / "repo/antarctica/results/timing"
    timing.mkdir(parents=True)
    mumps = "scpc_mumps_10step_dt0p125at2500_dg0_logvelnet_cached_strict_v4"
    (timing / "latest_matrix_campaign.txt").write_text(mumps + "\n")
    out = sandbox / "matrix.md"
    proc, _ = make(sandbox, "matrix", f"MATRIX_OUTPUT={out}")
    assert proc.returncode == 0, proc.stderr
    assert f"Campaign tag: `{mumps}`" in out.read_text()
    proc, _ = make(sandbox, "matrix", f"MATRIX_OUTPUT={out}", "TIMING_SOLVER=scpc_gamg")
    assert proc.returncode == 0, proc.stderr
    text = out.read_text()
    assert f"Campaign tag: `{mumps.replace('mumps', 'gamg')}`" in text
    assert "Lanes use `scpc_gamg`, an exact-mesh prepared cache (prepared under `scpc_mumps`" in text


def _python_shim(sandbox):
    import sys
    (sandbox / "bin" / "python").write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    (sandbox / "bin" / "python").chmod(0o755)


def _map_check_inputs(sandbox):
    r"""A fetched stand-in MAP, the production mesh and the two sidecars."""
    import hashlib
    import shutil
    root = sandbox / "repo" / "antarctica"
    basename = "inversion_icepack2_rc_n3_dg0_logvelnet_2000_int5000_bilap_snap20260922_0948.h5"
    download = root / "results/map_check/maps/download" / basename
    download.parent.mkdir(parents=True)
    download.write_text("a released MAP\n")
    (root / "mesh").mkdir()
    for name in ("boundary_ids_antarctica_10000_1000_buffered20000.json",
                 "boundary_ids_antarctica_5000_2000_buffered0.json"):
        shutil.copy(REPO / "antarctica" / "mesh" / name, root / "mesh" / name)
    (root / "mesh" / "antarctica_10000_1000_buffered20000.msh").write_text("mesh\n")
    return [f"MAP_CHECK_MD5={hashlib.md5(download.read_bytes()).hexdigest()}",
            "MAP_CHECK_ATTRS=0"]


def test_map_check_dry_run_lists_every_stage_and_submits_nothing(sandbox):
    _python_shim(sandbox)
    proc, seen = make(sandbox, "map-check-dry-run", *_map_check_inputs(sandbox),
                      ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    out = proc.stdout
    stages = ("fetch", "repack", "score_native", "prepare_transfer", "audit_cache",
              "score_transfer", "lane_transfer", "lane_native", "control_transfer",
              "control_native", "audit_controls", "summary")
    positions = [out.index(f"--- {stage}") for stage in stages]
    assert positions == sorted(positions)
    assert out.count("DRY RUN: site iu_quartz: sbatch") == 10
    assert "ISMIP7_FRICTION=regularized_coulomb" in out
    assert "ISMIP7_MESH=checkpoint" in out and "antarctica_10000_1000_buffered20000.msh" in out
    assert "--ntasks-per-node=16" in out and "--ntasks-per-node=64" in out
    assert "scripts/batch_runners/timing_prepare.script" in out
    assert "antarctica/scripts/batch_runners/projection.sbatch" in out
    assert seen == []


def test_map_check_submits_the_first_runnable_stage_and_stamps_it(sandbox):
    _python_shim(sandbox)
    inputs = _map_check_inputs(sandbox)
    # First call: the fetch runs in the manager itself and passes on the md5.
    proc, seen = make(sandbox, "map-check", *inputs, ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert "FETCH PASSED" in proc.stdout and seen == []
    # Second call: the repack is the one runnable stage and goes to sbatch.
    proc, seen = make(sandbox, "map-check", *inputs, ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert seen[-1] == "ARG: scripts/batch_runners/timing_redistribute.script"
    stem_dir = (sandbox / "repo/antarctica/results/map_check"
                / "inversion_icepack2_rc_n3_dg0_logvelnet_2000_int5000_bilap_snap20260922_0948")
    assert (stem_dir / "status_repack.txt").read_text().startswith("submitted job_id=5151 ")
    assert (stem_dir / "summary.md").is_file()
    # Third call: the repack is active and nothing else can run yet.
    (sandbox / "sbatch_calls.txt").unlink()
    proc, seen = make(sandbox, "map-check", *inputs, ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert seen == [] and "active: repack" in proc.stdout
