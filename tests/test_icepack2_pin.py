r"""`install_deps.sh --check-icepack2`: one known icepack2 tree at every site.

The inversion needs two lines Firedrake 2026 forces on icepack2, and they lived
as uncommitted edits in one workstation checkout, rsynced to the clusters, so no
reader could tell whether two sites ran the same code (issue #46). The source is
now a pinned commit, and these hold what the pin promises: a site with no
checkout gets the pinned one, a site already on it is left alone, and a site
carrying anything else is told rather than overwritten, since the rsynced copy
may be the only place its edits exist.

The pinned default is a real remote, so these run against a local repository
through the three ICEPACK2_* overrides instead of reaching the network. The
check needs no venv, but it goes through the site layer, so it is given a site
and a nominal ISMIP7_FIREDRAKE like the other runner tests.
"""
import os
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
INSTALL = REPO / "antarctica" / "scripts" / "batch_runners" / "install_deps.sh"


def git(*args, cwd):
    return subprocess.run(("git",) + args, cwd=cwd, check=True,
                          capture_output=True, text=True).stdout.strip()


@pytest.fixture
def origin(tmp_path):
    r"""A repository standing in for the pinned branch, and its head SHA."""
    d = tmp_path / "icepack2-origin"
    d.mkdir()
    git("init", "--quiet", "-b", "main", cwd=d)
    git("config", "user.email", "t@example.org", cwd=d)
    git("config", "user.name", "T", cwd=d)
    (d / "setup.py").write_text("# icepack2\n")
    git("add", "setup.py", cwd=d)
    git("commit", "--quiet", "-m", "base", cwd=d)
    base = git("rev-parse", "HEAD", cwd=d)
    git("checkout", "--quiet", "-b", "fix/shim", cwd=d)
    (d / "shim.py").write_text("d = mesh.geometric_dimension\n")
    git("add", "shim.py", cwd=d)
    git("commit", "--quiet", "-m", "shim", cwd=d)
    return d, base, git("rev-parse", "HEAD", cwd=d)


def check(work, origin_dir, sha, ref="fix/shim"):
    env = dict(os.environ, ISMIP7_SITE="local", ISMIP7_REPO=str(REPO),
               ISMIP7_WORK=str(work), FD_PREFIX=str(work / "sw"),
               ISMIP7_FIREDRAKE=str(INSTALL),
               ICEPACK2_REMOTE=str(origin_dir), ICEPACK2_REF=ref,
               ICEPACK2_SHA=sha)
    env.pop("VIRTUAL_ENV", None)
    return subprocess.run(["bash", str(INSTALL), "--check-icepack2"],
                          cwd=REPO, env=env, capture_output=True, text=True)


def test_a_site_with_no_checkout_gets_the_pinned_commit(tmp_path, origin):
    origin_dir, _, sha = origin
    proc = check(tmp_path, origin_dir, sha)
    assert proc.returncode == 0, proc.stderr
    assert sha in proc.stdout
    src = tmp_path / "sw" / "src" / "icepack2"
    assert git("rev-parse", "HEAD", cwd=src) == sha
    assert (src / "shim.py").exists()


def test_a_checkout_on_another_commit_is_moved_onto_the_pin(tmp_path, origin):
    origin_dir, base, sha = origin
    src = tmp_path / "sw" / "src" / "icepack2"
    src.parent.mkdir(parents=True)
    subprocess.run(["git", "clone", "--quiet", str(origin_dir), str(src)],
                   check=True)
    git("checkout", "--quiet", "--detach", base, cwd=src)
    proc = check(tmp_path, origin_dir, sha)
    assert proc.returncode == 0, proc.stderr
    assert git("rev-parse", "HEAD", cwd=src) == sha


def test_local_changes_are_reported_and_left_alone(tmp_path, origin):
    origin_dir, _, sha = origin
    src = tmp_path / "sw" / "src" / "icepack2"
    src.parent.mkdir(parents=True)
    subprocess.run(["git", "clone", "--quiet", str(origin_dir), str(src)],
                   check=True)
    (src / "setup.py").write_text("# edited by hand\n")
    proc = check(tmp_path, origin_dir, sha)
    assert proc.returncode == 2
    assert "local changes" in proc.stderr
    assert (src / "setup.py").read_text() == "# edited by hand\n"


def test_an_rsynced_copy_is_named_not_overwritten(tmp_path, origin):
    origin_dir, _, sha = origin
    src = tmp_path / "sw" / "src" / "icepack2"
    src.mkdir(parents=True)
    (src / "setup.py").write_text("# rsynced\n")
    proc = check(tmp_path, origin_dir, sha)
    assert proc.returncode == 2
    assert "not a git checkout" in proc.stderr and ".rsynced" in proc.stderr
    assert (src / "setup.py").read_text() == "# rsynced\n"


def test_a_commit_the_remote_does_not_carry_is_an_error(tmp_path, origin):
    origin_dir, _, _ = origin
    absent = "0" * 40
    proc = check(tmp_path, origin_dir, absent)
    assert proc.returncode == 2
    assert absent in proc.stderr or "fatal" in (proc.stderr + proc.stdout)


def test_the_default_pin_is_the_branch_of_the_open_pull_request():
    r"""The defaults are the record of which tree every site installs, so they
    are pinned here too: a change to them is a deliberate one, made when
    icepack/icepack2#3 merges and the fix can come from main."""
    text = INSTALL.read_text()
    assert 'ICEPACK2_REMOTE="${ICEPACK2_REMOTE:-https://github.com/hoffmaao/icepack2.git}"' in text
    assert 'ICEPACK2_REF="${ICEPACK2_REF:-fix/firedrake-2026-geometric-dimension}"' in text
    assert 'ICEPACK2_SHA="${ICEPACK2_SHA:-e0a46c9ce3e95dc916660a200e695f0417aa3037}"' in text
