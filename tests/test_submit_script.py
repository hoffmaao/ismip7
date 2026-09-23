r"""`submit.sh script`: the one place the timing campaign's sbatch line is made.

antarctica/Makefile and manage_timing_campaign.py both submit through it, so
these pin what they rely on: the site's account, node feature and extra flags
reach the command; standard output is sbatch's alone, because the callers read
the job id from it; the submission really happens from antarctica/, because the
timing scripts resolve scripts/, mesh/ and their logs from SLURM_SUBMIT_DIR;
the job's exit status comes back under --wait; and a request one node of the
site cannot hold is refused with status 3 and a single-token reason, which the
campaign records as not_runnable. Every KEY=VALUE travels in sbatch's own
environment under a bare --export=ALL, since sbatch splits an --export list on
commas and has Slurm rebuild the login environment at the start of a job given
one, and the printed line is the command that runs.
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
env | grep '^ISMIP7_' | sort | sed 's/^/ENV: /' >> "$SBATCH_CALLS"
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


def environment(seen):
    r"""The ISMIP7_* environment the fake sbatch ran with."""
    return dict(line[len("ENV: "):].split("=", 1) for line in seen if line.startswith("ENV: "))


def printed(output):
    r"""The composed command split as the shell splits it: the assignments
    `env` makes, then sbatch's own argv."""
    words = shlex.split(output.strip().split(": ", 1)[1])
    assert words[0] == "env", words
    return words[1:words.index("sbatch")], words[words.index("sbatch"):]


def test_quartz_line_carries_the_account_and_the_lane_request(bin_dir):
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--queue", "short",
                  "--tasks", "32", "--mem", "240G", "--time", "12:00:00",
                  "--name", "timing_500_5000_32", "--dependency", "afterok:7",
                  "--dry-run", "ISMIP7_LC=500", ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""
    assert proc.stderr.startswith("site iu_quartz: env ISMIP7_SITE=iu_quartz ")
    assignments, argv = printed(proc.stderr)
    assert assignments == ["ISMIP7_SITE=iu_quartz", f"ISMIP7_REPO={REPO}", "ISMIP7_LC=500"]
    assert argv[:8] == ["sbatch", "-A", "r00905", "--parsable",
                        "-J", "timing_500_5000_32", "-p", "general"]
    for flag in ("--nodes=1", "--ntasks-per-node=32", "--mem=240G", "--time=12:00:00",
                 "--dependency=afterok:7", "--export=ALL"):
        assert flag in argv, argv
    assert argv[-1] == SCRIPT
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
    assert proc.stderr.startswith("site local: env ISMIP7_SITE=local ")
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


def test_every_value_reaches_sbatch_s_environment_under_a_bare_all(bin_dir):
    r"""An --export list has two faults. sbatch splits it on commas, so
    ISMIP7_SUBCYCLES=1,4,16,64 would reach the job as ISMIP7_SUBCYCLES=1, and
    any list sets SLURM_GET_USER_ENV=1, under which Slurm rebuilds the login
    environment when the job starts and holds the job when that fails. In
    sbatch's own environment a value arrives whole, an empty one stays set,
    and each wins over the same variable exported in the calling shell."""
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica",
                  "ISMIP7_SUBCYCLES=1,4,16,64", "ISMIP7_LC=500", "ISMIP7_APPARENT_MB=",
                  ISMIP7_SITE="local", ISMIP7_FIREDRAKE=str(SUBMIT), ISMIP7_SUBCYCLES="1")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == "4242;cluster\n"
    seen = calls(bin_dir)
    env = environment(seen)
    assert env["ISMIP7_SUBCYCLES"] == "1,4,16,64", seen
    assert env["ISMIP7_LC"] == "500" and env["ISMIP7_APPARENT_MB"] == ""
    assert [line for line in seen if line.startswith("ARG: --export")] == ["ARG: --export=ALL"]
    # env execs sbatch, so under --wait the job's status still comes back.
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica", "--wait",
                  "ISMIP7_SUBCYCLES=1,4,16,64", ISMIP7_SITE="local",
                  ISMIP7_FIREDRAKE=str(SUBMIT), FAKE_SBATCH_RC="7")
    assert proc.returncode == 7


def test_the_printed_line_is_the_command_that_runs(bin_dir):
    r"""Split as the shell splits it, the printed line gives the environment
    sbatch ran with and the arguments it received, the site and the checkout
    first."""
    proc = submit(bin_dir, "script", SCRIPT, "--cd", "antarctica",
                  "ISMIP7_SUBCYCLES=1,4,16,64", "ISMIP7_EXPERIMENT=hist_cesm_waccm",
                  ISMIP7_SITE="local", ISMIP7_FIREDRAKE=str(SUBMIT))
    assert proc.returncode == 0, proc.stderr
    assignments, argv = printed(proc.stderr)
    assert assignments == ["ISMIP7_SITE=local", f"ISMIP7_REPO={REPO}",
                           "ISMIP7_SUBCYCLES=1,4,16,64", "ISMIP7_EXPERIMENT=hist_cesm_waccm"]
    seen = calls(bin_dir)
    assert [line[len("ARG: "):] for line in seen if line.startswith("ARG: ")] == argv[1:]
    env = environment(seen)
    for name, value in (assignment.split("=", 1) for assignment in assignments):
        assert env[name] == value, (name, seen)
    # The runners' own kinds print the same form, on standard output.
    proc = submit(bin_dir, "projection", "--dry-run", "ISMIP7_SUBCYCLES=1,4,16,64",
                  ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
    assignments, argv = printed(proc.stdout)
    assert assignments == ["ISMIP7_SITE=iu_quartz", f"ISMIP7_REPO={REPO}",
                           "ISMIP7_SUBCYCLES=1,4,16,64"]
    assert [word for word in argv if word.startswith("--export")] == ["--export=ALL"]
    assert argv[-1] == "antarctica/scripts/batch_runners/projection.sbatch"


def test_a_key_that_is_no_variable_name_is_refused(bin_dir):
    proc = submit(bin_dir, "projection", "--dry-run", "ISMIP7-LC=500",
                  ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 2 and "ISMIP7-LC" in proc.stderr
    assert proc.stdout == ""


def test_several_pairs_in_one_argument_are_refused(bin_dir):
    r"""zsh passes an unquoted $VAR as one argument, so `submit.sh projection
    $COMMON` would set ISMIP7_LC to "32000 ISMIP7_LC_COARSE=320000 ..."."""
    joined = "ISMIP7_LC=32000 ISMIP7_LC_COARSE=320000 ISMIP7_SUBCYCLES=1,4,16,64"
    proc = submit(bin_dir, "projection", "--dry-run", joined, ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 2 and "several KEY=VALUE pairs" in proc.stderr
    assert "${=VAR}" in proc.stderr
    assert proc.stdout == ""
    # A space alone, or an option after one, still makes a single value.
    proc = submit(bin_dir, "projection", "--dry-run", "ISMIP7_RUN_TAG=a b",
                  "ISMIP7_NOTE=--mail-type=FAIL --mail-user=me@example.org",
                  ISMIP7_SITE="iu_quartz")
    assert proc.returncode == 0, proc.stderr
