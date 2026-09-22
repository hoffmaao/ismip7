r"""How the timing campaign hands a lane to the scheduler.

`CampaignManager._submit` goes through `batch_runners/submit.sh script`, so a
lane is the same request at every site and the site file supplies the account,
node feature and per-node limits. These drive the real submit.sh against a fake
sbatch and check the three things the campaign's state machine depends on: the
job id it stamps, the `not_runnable` record for a lane this site's nodes cannot
hold, and a dry run that works where no site matches.
"""
import os
import sys
from pathlib import Path

import pytest

pytest.importorskip("firedrake")

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "antarctica"
sys.path.insert(0, str(ROOT / "scripts"))

import manage_timing_campaign as mtc  # noqa: E402

SCRIPT = ROOT / "scripts/batch_runners/timing_transient.script"

FAKE_SBATCH = '''#!/bin/bash
printf 'CWD: %s\\n' "$PWD" >> "$SBATCH_CALLS"
env | grep '^ISMIP7_' | sort | sed 's/^/ENV: /' >> "$SBATCH_CALLS"
printf 'ARG: %s\\n' "$@" >> "$SBATCH_CALLS"
[ "${FAKE_SBATCH_RC:-0}" = 0 ] || { echo "sbatch: error: Batch job submission failed" >&2; exit "$FAKE_SBATCH_RC"; }
echo "4242;cluster"
'''


@pytest.fixture
def manager(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "sbatch").write_text(FAKE_SBATCH)
    (bin_dir / "hostname").write_text("#!/bin/bash\necho nowhere.example\n")
    for stub in bin_dir.iterdir():
        stub.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{os.environ['PATH']}")
    monkeypatch.setenv("ISMIP7_LOCAL_ENV", os.devnull)
    monkeypatch.setenv("SBATCH_CALLS", str(tmp_path / "sbatch_calls.txt"))
    monkeypatch.setenv("ISMIP7_FIREDRAKE", str(SCRIPT))   # any readable file
    for name in ("ISMIP7_SITE", "ISMIP7_REPO", "ISMIP7_CORES_PER_NODE",
                 "ISMIP7_MEM_PER_NODE", "FAKE_SBATCH_RC"):
        monkeypatch.delenv(name, raising=False)
    # _submit reads only these; the constructor pins campaign parameters into
    # os.environ, which a unit test has no business doing.
    m = mtc.CampaignManager.__new__(mtc.CampaignManager)
    m.root = ROOT
    m.queue, m.partition, m.constraint = "short", None, None
    m.walltime = "12:00:00"
    m.dry_run = False
    m.submit_failures = 0
    m.tmp = tmp_path
    return m


def submit(m, ncores=16, memory="64G", **kwargs):
    status = m.tmp / "status.txt"
    m._submit("timing_2500_25000_16", ncores, memory,
              {"ISMIP7_LC": 2500, "ISMIP7_TIMING_STATUS": status},
              SCRIPT, status, **kwargs)
    return mtc.read_status(status)


def calls(m):
    path = m.tmp / "sbatch_calls.txt"
    return path.read_text().splitlines() if path.exists() else []


def test_a_lane_is_submitted_from_antarctica_and_its_job_id_stamped(manager, monkeypatch, capsys):
    monkeypatch.setenv("ISMIP7_SITE", "iu_quartz")
    status = submit(manager, dependency="afterok:7")
    assert status["state"] == "submitted" and status["job_id"] == "4242"
    seen = calls(manager)
    assert seen[0] == f"CWD: {ROOT}"
    assert seen[-1] == "ARG: scripts/batch_runners/timing_transient.script"
    for arg in ("-A", "r00905", "-p", "general", "--ntasks-per-node=16", "--mem=64G",
                "--time=12:00:00", "--dependency=afterok:7", "timing_2500_25000_16"):
        assert f"ARG: {arg}" in seen, seen
    # Every value travels in sbatch's own environment under a bare ALL.
    assert "ARG: --export=ALL" in seen
    assert "ENV: ISMIP7_LC=2500" in seen and f"ENV: ISMIP7_REPO={REPO}" in seen
    assert "SUBMIT:" in capsys.readouterr().out


def test_a_partition_by_name_replaces_the_queue_class(manager, monkeypatch):
    monkeypatch.setenv("ISMIP7_SITE", "iu_quartz")
    manager.partition = "debug"
    submit(manager)
    assert "ARG: debug" in calls(manager) and "ARG: general" not in calls(manager)


def test_a_lane_the_site_s_nodes_cannot_hold_is_recorded_not_runnable(manager, monkeypatch, capsys):
    monkeypatch.setenv("ISMIP7_SITE", "rice_nots")
    monkeypatch.setenv("ISMIP7_CORES_PER_NODE", "40")
    monkeypatch.setenv("ISMIP7_MEM_PER_NODE", "187G")
    status = submit(manager, ncores=64)
    assert status["state"] == "not_runnable"
    assert status["reason"] == "exceeds_site_cores"
    assert calls(manager) == [] and manager.submit_failures == 0
    status = submit(manager, ncores=32, memory="240G")
    assert status["reason"] == "exceeds_site_mem"
    assert "NOT RUNNABLE:" in capsys.readouterr().out
    # A node feature named by hand is the operator's claim that it fits.
    manager.constraint = "sapphirerapids"
    assert submit(manager, ncores=64, memory="240G")["state"] == "submitted"
    assert "ARG: sapphirerapids" in calls(manager)


def test_a_refused_submission_is_a_failure(manager, monkeypatch):
    monkeypatch.setenv("ISMIP7_SITE", "iu_quartz")
    monkeypatch.setenv("FAKE_SBATCH_RC", "1")
    status = submit(manager)
    assert status["state"] == "submission_failed" and status["exit_code"] == "1"
    assert manager.submit_failures == 1


def test_a_dry_run_works_where_no_site_matches_and_stamps_nothing(manager, capsys):
    manager.dry_run = True
    assert submit(manager) is None
    out = capsys.readouterr().out
    assert out.startswith("DRY RUN: site local: env ISMIP7_SITE=local ")
    assert "--ntasks-per-node=16" in out
    assert calls(manager) == [] and manager.submit_failures == 0


def test_a_real_submission_never_guesses_a_site(manager):
    status = submit(manager)
    assert status["state"] == "submission_failed"
    assert calls(manager) == []


def test_a_record_says_where_it_was_measured_and_from_what_cache(tmp_path):
    r"""seconds_per_step has no warm-up excluded, and the matrix is run at
    several sites. The record carries the site and whether the kernel cache
    the lane started from was empty; none of it reaches a cache manifest."""
    import timing_campaign as tc
    cache = tmp_path / "pyop2"
    cache.mkdir()
    cold = tc.host_provenance({"ISMIP7_SITE": "rice_nots", "SLURM_JOB_ID": "12",
                               "PYOP2_CACHE_DIR": str(cache)})
    assert cold["site"] == "rice_nots" and cold["slurm_job_id"] == "12"
    assert cold["jit_cache_dir"] == str(cache) and cold["jit_cache_was_empty"] is True
    (cache / "kernel.so").write_text("x")
    assert tc.host_provenance({"PYOP2_CACHE_DIR": str(cache)})["jit_cache_was_empty"] is False
    default = tc.host_provenance({})
    assert default["jit_cache_dir"] is None and default["jit_cache_was_empty"] is None
    assert default["site"] is None and default["container"] is None
    assert tc.CAMPAIGN_VERSION == 4
