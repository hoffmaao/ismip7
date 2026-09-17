r"""The site half of the batch environment, exercised without a scheduler.

`site_core.sh` picks a cluster definition, reads the per-user `sites/local.env`
ahead of it, and must not choose a model configuration: the timing campaign's
job scripts source it alone, `run_timing.py` refuses a matrix lane that arrives
with `ISMIP7_MESH` set, and a defaulted `ISMIP7_FRICTION` would change the law
under any stage that did not export one. `site_env.sh` adds those defaults for
the inversion and projection runners.

Each test sources one of the two files in a clean bash and reads variables back.
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BR = REPO / "antarctica" / "scripts" / "batch_runners"

MODEL_VARS = (
    "ISMIP7_FRICTION", "ISMIP7_MESH", "ISMIP7_LC", "ISMIP7_LC_COARSE",
    "ISMIP7_GEOMETRY_SPACE", "ISMIP7_N_FLOW", "ISMIP7_MAP_DEFAULT",
    "ISMIP7_DATA_ROOT", "ISMIP7_OBS_DATA_ROOT",
)

# What site_env.sh hangs off ISMIP7_SHARE: the artifacts a checkout does not
# carry, as opposed to ISMIP7_REPO's code.
SHARED_VARS = (
    "ISMIP7_DATA_ROOT", "ISMIP7_OBS_DATA_ROOT", "ISMIP7_MESH", "ISMIP7_MAP_DEFAULT",
)


@pytest.fixture
def runners(tmp_path):
    r"""A private copy of batch_runners, so a test can write sites/local.env
    and add site files without touching the checkout."""
    dest = tmp_path / "repo" / "antarctica" / "scripts" / "batch_runners"
    shutil.copytree(BR, dest, ignore=shutil.ignore_patterns("local.env"))
    shutil.copy(BR.parent / "ismip7_names.sh", dest.parent / "ismip7_names.sh")
    return dest


def source(runners, file, expr, **env):
    r"""Source `file` and echo `expr`. Returns (exit code, stdout, stderr)."""
    base = {"PATH": os.environ["PATH"], "HOME": str(runners)}
    base.update(env)
    proc = subprocess.run(
        ["bash", "-c", f'. "{runners / file}"; {expr}'],
        env=base, capture_output=True, text=True,
    )
    return proc.returncode, proc.stdout.strip(), proc.stderr


def show(*names):
    return "; ".join(f'echo "{n}=${{{n}-<unset>}}"' for n in names)


@pytest.mark.parametrize("site", ["local", "iu_quartz", "rice_nots", "uchicago_midway"])
def test_the_core_file_chooses_no_model_configuration(runners, site):
    r"""Every site, not just the no-scheduler one: a site file that named a mesh
    would defeat the timing lanes' `ISMIP7_MESH is required` guard, since they
    source this file alone and bring their own model configuration."""
    rc, out, err = source(runners, "site_core.sh", show(*MODEL_VARS),
                          ISMIP7_SITE=site)
    assert rc == 0, err
    assert out.splitlines() == [f"{n}=<unset>" for n in MODEL_VARS]


def test_site_env_still_supplies_the_runners_model_defaults(runners):
    rc, out, err = source(runners, "site_env.sh", show(*MODEL_VARS),
                          ISMIP7_SITE="local", ISMIP7_REPO="/repo")
    assert rc == 0, err
    values = dict(line.split("=", 1) for line in out.splitlines())
    assert values["ISMIP7_FRICTION"] == "regularized_coulomb"
    assert values["ISMIP7_LC"] == "2500"
    assert values["ISMIP7_MESH"] == "/repo/antarctica/mesh/antarctica_64000_2500.msh"
    assert values["ISMIP7_DATA_ROOT"] == "/repo/ISMIP7/AIS"
    assert values["ISMIP7_OBS_DATA_ROOT"] == "/repo/antarctica/data"
    assert "<unset>" not in out


def test_the_banner_prints_the_model_lines_only_with_site_env(runners):
    env = dict(ISMIP7_SITE="local", ISMIP7_REPO="/repo")
    rc, core, err = source(runners, "site_core.sh", "ismip7_banner", **env)
    assert rc == 0, err
    assert "repo    /repo" in core
    assert "friction=" not in core and "mesh " not in core
    rc, full, err = source(runners, "site_env.sh", "ismip7_banner", **env)
    assert rc == 0, err
    assert full.splitlines()[-2:] == [
        "    mesh    /repo/antarctica/mesh/antarctica_64000_2500.msh",
        "    lc=2500 lc_coarse=64000 geometry=dg0 friction=regularized_coulomb n=3.0",
    ]


@pytest.mark.parametrize("site", ["iu_quartz", "rice_nots", "uchicago_midway"])
def test_repo_self_is_the_checkout_the_file_was_sourced_from(runners, site):
    r"""Every cluster site, not just one. submit.sh cds to ISMIP7_REPO, so a
    site that pins it submits that tree's code from any other checkout -- which
    is how rice_nots kept pointing at one September checkout while this test
    covered only Quartz."""
    rc, out, err = source(runners, "site_core.sh", show("ISMIP7_REPO_SELF", "ISMIP7_REPO"),
                          ISMIP7_SITE=site)
    assert rc == 0, err
    repo = str(runners.parents[2])
    assert out.splitlines() == [f"ISMIP7_REPO_SELF={repo}", f"ISMIP7_REPO={repo}"]


def test_a_site_that_states_no_share_keeps_the_artifacts_beside_its_code(runners):
    r"""ISMIP7_SHARE defaults to ISMIP7_REPO, so a cluster holding one checkout
    is unaffected by the split."""
    rc, out, err = source(runners, "site_env.sh", show("ISMIP7_SHARE", *SHARED_VARS),
                          ISMIP7_SITE="local", ISMIP7_REPO="/repo")
    assert rc == 0, err
    values = dict(line.split("=", 1) for line in out.splitlines())
    assert values["ISMIP7_SHARE"] == "/repo"
    assert all(values[n].startswith("/repo/") for n in SHARED_VARS), values


def test_rice_keeps_the_shared_artifacts_off_the_checkout(runners):
    r"""The converse of the test above: the code root follows the invocation,
    but the ~313 GB forcing tree, the shared meshes, the MAPs and the
    observational rasters do not exist in a second checkout, so they must stay
    on the tree that holds them."""
    rc, out, err = source(runners, "site_env.sh",
                          show("ISMIP7_REPO", "ISMIP7_SHARE", *SHARED_VARS),
                          ISMIP7_SITE="rice_nots")
    assert rc == 0, err
    values = dict(line.split("=", 1) for line in out.splitlines())
    checkout = str(runners.parents[2])
    assert values["ISMIP7_REPO"] == checkout
    assert values["ISMIP7_SHARE"] == "/projects/ah301/ismip7"
    for name in SHARED_VARS:
        assert not values[name].startswith(checkout), name
        assert values[name].startswith("/projects/ah301/ismip7/"), name


def test_a_share_given_for_one_submission_moves_every_artifact(runners):
    r"""One variable answers for all four, so a second share needs one word."""
    rc, out, err = source(runners, "site_env.sh", show(*SHARED_VARS),
                          ISMIP7_SITE="rice_nots", ISMIP7_SHARE="/tmp/share")
    assert rc == 0, err
    values = dict(line.split("=", 1) for line in out.splitlines())
    assert all(values[n].startswith("/tmp/share/") for n in SHARED_VARS), values


@pytest.mark.parametrize("name", SHARED_VARS)
def test_a_shared_path_given_for_one_submission_still_wins(runners, name):
    rc, out, err = source(runners, "site_env.sh", show(name),
                          ISMIP7_SITE="rice_nots", **{name: "/tmp/elsewhere"})
    assert rc == 0, err
    assert out == f"{name}=/tmp/elsewhere"


def test_repo_self_stays_inside_a_symlinked_sandbox(runners, tmp_path):
    r"""The chain tests reach batch_runners through a symlink. Resolving it
    would point ISMIP7_REPO at the real checkout rather than the sandbox."""
    sandbox = tmp_path / "sandbox"
    sandbox.mkdir()
    (sandbox / "antarctica").symlink_to(runners.parents[1])
    linked = sandbox / "antarctica" / "scripts" / "batch_runners"
    rc, out, err = source(linked, "site_core.sh", show("ISMIP7_REPO_SELF"),
                          ISMIP7_SITE="local")
    assert rc == 0, err
    assert out == f"ISMIP7_REPO_SELF={sandbox}"


def test_local_env_beats_the_site_file_and_loses_to_the_command_line(runners):
    (runners / "sites" / "local.env").write_text(
        'ISMIP7_SITE="${ISMIP7_SITE:-iu_quartz}"\n'
        'ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-mine}"\n'
    )
    expr = show("ISMIP7_SITE_NAME", "ISMIP7_ACCOUNT", "ISMIP7_PART_DEBUG")
    rc, out, err = source(runners, "site_core.sh", expr)
    assert rc == 0, err
    assert out.splitlines() == [
        "ISMIP7_SITE_NAME=iu_quartz", "ISMIP7_ACCOUNT=mine", "ISMIP7_PART_DEBUG=debug",
    ]
    rc, out, err = source(runners, "site_core.sh", expr,
                          ISMIP7_ACCOUNT="once", ISMIP7_SITE="rice_nots")
    assert rc == 0, err
    assert out.splitlines() == [
        "ISMIP7_SITE_NAME=rice_nots", "ISMIP7_ACCOUNT=once", "ISMIP7_PART_DEBUG=scavenge",
    ]


def test_local_env_override_names_another_file(runners, tmp_path):
    (runners / "sites" / "local.env").write_text('ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-mine}"\n')
    rc, out, err = source(runners, "site_core.sh", show("ISMIP7_ACCOUNT"),
                          ISMIP7_SITE="rice_nots", ISMIP7_LOCAL_ENV=os.devnull)
    assert rc == 0, err
    assert out == "ISMIP7_ACCOUNT="


def test_local_env_is_never_taken_for_a_site_definition(runners, tmp_path):
    r"""The hostname loop reads ISMIP7_SITE_MATCH out of sites/*.sh. A local.env
    that mentions one (copied from a site file, say) must not become a site."""
    (runners / "sites" / "local.env").write_text('ISMIP7_SITE_MATCH="*"\n')
    # A host no shipped site claims, wherever the test itself runs.
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "hostname").write_text("#!/bin/bash\necho nowhere.example\n")
    (bin_dir / "hostname").chmod(0o755)
    rc, out, err = source(runners, "site_core.sh", "true",
                          PATH=f"{bin_dir}:{os.environ['PATH']}")
    assert rc == 2
    assert "no site definition matches host" in err


def test_a_site_with_an_empty_required_field_names_it(runners):
    rc, out, err = source(runners, "site_core.sh", "ismip7_site_require",
                          ISMIP7_SITE="uchicago_midway")
    assert rc == 2
    assert "is missing:" in err and "ISMIP7_PART_LONG" in err and "ISMIP7_WORK" in err
    assert "ISMIP7_REPO" not in err.split("is missing:")[1].splitlines()[0]


@pytest.mark.parametrize("env,expected", [
    (dict(ISMIP7_SITE="iu_quartz"), ["", "128", "515700M"]),
    (dict(ISMIP7_SITE="rice_nots"), ["cascadelake", "40", "187G"]),
    (dict(ISMIP7_SITE="rice_nots", ISMIP7_CONSTRAINT_TIMING="sapphirerapids"),
     ["sapphirerapids", "96", ""]),
    (dict(ISMIP7_SITE="local"), ["", "", ""]),
])
def test_the_per_node_limits_describe_the_timing_node_class(runners, env, expected):
    r"""`submit.sh script` refuses what one node cannot hold, so the limits have
    to be those of the nodes ISMIP7_CONSTRAINT_TIMING selects, and a site that
    states none (the no-scheduler one) refuses nothing."""
    names = ("ISMIP7_CONSTRAINT_TIMING", "ISMIP7_CORES_PER_NODE", "ISMIP7_MEM_PER_NODE")
    rc, out, err = source(runners, "site_core.sh", show(*names), **env)
    assert rc == 0, err
    assert [line.split("=", 1)[1] for line in out.splitlines()] == expected


@pytest.mark.parametrize("host,site", [
    ("h1.quartz.uits.iu.edu", "iu_quartz"),     # the login node the pattern used to miss
    ("h2.quartz.uits.iu.edu", "iu_quartz"),
    ("c42.quartz.uits.iu.edu", "iu_quartz"),
    ("login3.nots.rice.edu", "rice_nots"),
    ("midway3-login4.rcc.uchicago.edu", "uchicago_midway"),
])
def test_each_cluster_s_own_hostnames_choose_its_site(runners, tmp_path, host, site):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "hostname").write_text(f"#!/bin/bash\necho {host}\n")
    (bin_dir / "hostname").chmod(0o755)
    rc, out, err = source(runners, "site_core.sh", show("ISMIP7_SITE_NAME"),
                          PATH=f"{bin_dir}:{os.environ['PATH']}")
    assert rc == 0, err
    assert out == f"ISMIP7_SITE_NAME={site}"
