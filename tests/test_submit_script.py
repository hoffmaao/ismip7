r"""`submit.sh script`: the one place the timing campaign's sbatch line is made.

antarctica/Makefile and manage_timing_campaign.py both submit through it, so
these pin what they rely on: the site's account, node feature and extra flags
reach the command; standard output is sbatch's alone, because the callers read
the job id from it; the submission really happens from antarctica/, because the
timing scripts resolve scripts/, mesh/ and their logs from SLURM_SUBMIT_DIR;
the job's exit status comes back under --wait; and a request one node of the
site cannot hold is refused with status 3 and a single-token reason, which the
campaign records as not_runnable. A KEY=VALUE whose value holds a comma travels
in sbatch's environment, since sbatch splits its --export list on commas, and so
reaches the job whole.
"""
import os
import shlex
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SUBMIT = REPO / "antarctica" / "scripts" / "batch_runners" / "submit.sh"
SCRIPT = "scripts/batch_runners/timing_transient.script"

FAKE_SBATCH = '''#!/bin/bash
printf 'CWD: %s\\n' "$PWD" >> "$SBATCH_CALLS"
[ -z "${ISMIP7_SUBCYCLES+x}" ] || printf 'ENV: ISMIP7_SUBCYCLES=%s\\n' "$ISMIP7_SUBCYCLES" >> "$SBATCH_CALLS"
printf 'ARG: %s\\n' "$@" >> "$SBATCH_CALLS"
echo "4242;cluster"
exit "${FAKE_SBATCH_RC:-0}"
'''


@pytest.fixture
def bin_dir(tmp_path):
    d = tmp_path / "bin"
    d.mkdir()
    (d / "sbatch").write_text(FAKE_SBATCH)
    (d / "sbatch").chmod(0o755)
    # A host no shipped site claims, wherever the test itself runs.
    (d / "hostname").write_text("#!/bin/bash\necho nowhere.example\n")
    (d / "hostname").chmod(0o755)
    return d


def submit(bin_dir, *args, **env):
    base = {
        "PATH": f"{bin_dir}:{os.environ['PATH']}",
        "HOME": str(bin_dir.parent),
        "ISMIP7_LOCAL_ENV": os.devnull,
        "ISMIP7_REPO": str(REPO),
        "SBATCH_CALLS": str(bin_dir.parent / "sbatch_calls.txt"),
    }
    base.update(env)
    return subprocess.run(["bash", str(SUBMIT), *args], env=base,
                          capture_output=True, text=True)


def calls(bin_dir):
    path = bin_dir.parent / "sbatch_calls.txt"
    return path.read_text().splitlines() if path.exists() else []


def test_quartz_line_carries_the_account_and_the_lane_request(bin_dir):
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--queue", "short",
                  "--tasks", "32", "--mem", "240G", "--time", "12:00:00",
                  "--name", "timing_500_5000_32", "--dependency", "afterok:7",
                  "--dry-run", "ISMIP7_LC=500", ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""
    line = proc.stderr.strip()
    assert line.startswith("site iu_quartz: sbatch -A r00905 --parsable -J timing_500_5000_32 -p general ")
    for flag in ("--nodes=1", "--ntasks-per-node=32", "--mem=240G", "--time=12:00:00",
                 "--dependency=afterok:7", "ISMIP7_LC=500", "ISMIP7_SITE=iu_quartz"):
        assert flag in line, line
    assert line.endswith(SCRIPT)
    assert calls(bin_dir) == []


def test_queue_classes_map_to_the_site_s_partitions(bin_dir):
    for queue, part in (("short", "commons"), ("long", "long"), ("debug", "scavenge")):
        proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica",
                      "--queue", queue, "--dry-run", ISMIP7_SITE="rice_nots")
        assert proc.returncode == 0, proc.stderr
        assert f" -p {part} " in proc.stderr
        # The timing node feature defaults to the forward's, and Rice needs no -A.
        assert "-C cascadelake" in proc.stderr and " -A " not in proc.stderr
    proc = submit(bin_dir, "script", SCRIPT, "--queue", "fast", "--dry-run",
                  ISMIP7_SITE="rice_nots")
    assert proc.returncode == 2


def test_site_extra_flags_precede_the_script(bin_dir):
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--wait", "--dry-run",
                  ISMIP7_SITE="iu_quartz",
                  ISMIP7_SBATCH_EXTRA="--mail-type=FAIL --mail-user=me@example.org")
    assert proc.returncode == 0, proc.stderr
    assert proc.stderr.strip().endswith(
        f"--mail-type=FAIL --mail-user=me@example.org --wait {SCRIPT}")
    # The runners' kinds get the same flags.
    proc = submit(bin_dir, "inversion", "--dry-run", ISMIP7_SITE="iu_quartz",
                  ISMIP7_SBATCH_EXTRA="--qos=long")
    assert proc.stdout.strip().endswith("--qos=long antarctica/scripts/batch_runners/inversion.sbatch")


def test_stdout_is_sbatch_s_and_the_job_starts_in_the_named_directory(bin_dir):
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--tasks", "16",
                  ISMIP7_SITE="local", ISMIP7_FIREDRAKE=str(SUBMIT))
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == "4242;cluster\n"
    assert "site local: sbatch" in proc.stderr
    seen = calls(bin_dir)
    assert seen[0] == f"CWD: {REPO / 'antarctica'}"
    assert seen[-1] == f"ARG: {SCRIPT}"
    assert not (REPO / "antarctica" / "logs").exists()


def test_a_failed_job_under_wait_fails_the_submission(bin_dir):
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--wait",
                  ISMIP7_SITE="local", ISMIP7_FIREDRAKE=str(SUBMIT), FAKE_SBATCH_RC="7")
    assert proc.returncode == 7
    assert "ARG: --wait" in calls(bin_dir)


@pytest.mark.parametrize("tasks,mem,reason", [
    ("64", "96G", "exceeds_site_cores"),
    ("32", "240G", "exceeds_site_mem"),
    ("32", "191489M", "exceeds_site_mem"),
])
def test_a_request_one_node_cannot_hold_is_refused_with_status_3(bin_dir, tasks, mem, reason):
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica",
                  "--tasks", tasks, "--mem", mem,
                  ISMIP7_SITE="rice_nots", ISMIP7_FIREDRAKE=str(SUBMIT),
                  ISMIP7_CORES_PER_NODE="40", ISMIP7_MEM_PER_NODE="187G")
    assert proc.returncode == 3
    assert f" reason={reason} " in proc.stderr
    assert proc.stdout == "" and calls(bin_dir) == []


def test_limits_admit_what_fits_and_yield_to_a_named_constraint(bin_dir):
    caps = dict(ISMIP7_SITE="rice_nots", ISMIP7_CORES_PER_NODE="40",
                ISMIP7_MEM_PER_NODE="187G")
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica",
                  "--tasks", "40", "--mem", "187G", "--dry-run", **caps)
    assert proc.returncode == 0, proc.stderr
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--tasks", "64",
                  "--mem", "240G", "--constraint", "sapphirerapids", "--dry-run", **caps)
    assert proc.returncode == 0, proc.stderr
    assert "-C sapphirerapids" in proc.stderr
    # The runners' own kinds are sized by the site file and are not checked.
    proc = submit(bin_dir, "inversion", "--dry-run", **caps)
    assert proc.returncode == 0, proc.stderr


def test_a_dry_run_needs_no_firedrake_but_a_submission_does(bin_dir):
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--dry-run",
                  ISMIP7_SITE="local")
    assert proc.returncode == 0, proc.stderr
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", ISMIP7_SITE="local")
    assert proc.returncode == 2 and "ISMIP7_FIREDRAKE" in proc.stderr
    assert calls(bin_dir) == []


def test_an_unknown_host_is_an_error_not_a_guess(bin_dir):
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--dry-run")
    assert proc.returncode == 2
    assert "no site definition matches host 'nowhere.example'" in proc.stderr


def test_script_options_are_refused_on_the_runners_kinds(bin_dir):
    proc = submit(bin_dir, "inversion", "--wait", "--dry-run", ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 2 and "submit.sh script" in proc.stderr


def test_a_missing_script_is_named(bin_dir):
    proc = submit(bin_dir, "script", "scripts/batch_runners/nope.script",
                  "--cd", "antarctica", "--dry-run", ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 2
    assert "antarctica/scripts/batch_runners/nope.script does not exist" in proc.stderr


def test_a_value_holding_a_comma_reaches_the_job_whole(bin_dir):
    r"""sbatch splits --export on commas, so in that list
    ISMIP7_SUBCYCLES=1,4,16,64 would reach the job as ISMIP7_SUBCYCLES=1.
    Set in sbatch's own environment, which ALL carries, it arrives whole, and
    it wins over the same variable exported in the calling shell."""
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica",
                  "ISMIP7_SUBCYCLES=1,4,16,64", "ISMIP7_LC=500",
                  ISMIP7_SITE="local", ISMIP7_FIREDRAKE=str(SUBMIT), ISMIP7_SUBCYCLES="1")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == "4242;cluster\n"
    seen = calls(bin_dir)
    assert "ENV: ISMIP7_SUBCYCLES=1,4,16,64" in seen, seen
    export = next(line for line in seen if line.startswith("ARG: --export="))
    fields = export.split("=", 1)[1].split(",")
    assert "ISMIP7_LC=500" in fields
    # Nothing in the list overrides the value or is a piece cut from it.
    assert not any(field.startswith("ISMIP7_SUBCYCLES=") for field in fields), fields
    assert all("=" in field for field in fields[1:]), fields
    # env execs sbatch, so under --wait the job's status still comes back.
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--wait",
                  "ISMIP7_SUBCYCLES=1,4,16,64", ISMIP7_SITE="local",
                  ISMIP7_FIREDRAKE=str(SUBMIT), FAKE_SBATCH_RC="7")
    assert proc.returncode == 7


def test_the_printed_line_sets_a_comma_value_ahead_of_sbatch(bin_dir):
    proc = submit(bin_dir, "projection", "--dry-run", "ISMIP7_SUBCYCLES=1,4,16,64",
                  "ISMIP7_EXPERIMENT=hist_cesm_waccm", ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    words = shlex.split(proc.stdout.split(": ", 1)[1])
    assert words[:3] == ["env", "ISMIP7_SUBCYCLES=1,4,16,64", "sbatch"]
    export = next(word for word in words if word.startswith("--export="))
    assert "ISMIP7_EXPERIMENT=hist_cesm_waccm" in export.split(",")
    assert "ISMIP7_SUBCYCLES" not in export
    assert words[-1] == "antarctica/scripts/batch_runners/projection.sbatch"
    assert calls(bin_dir) == []


def test_a_key_that_is_no_variable_name_is_refused(bin_dir):
    proc = submit(bin_dir, "projection", "--dry-run", "ISMIP7-LC=500",
                  ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 2 and "ISMIP7-LC" in proc.stderr
    assert proc.stdout == ""
