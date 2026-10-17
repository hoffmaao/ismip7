r"""ISMIP7_GRAD_CHECK=1 runs its Taylor test on the TAO path only.

On the default scipy metric (ISMIP7_GRAD_PRECOND=none) the knob used to be
ignored: the driver ran an ordinary L-BFGS-B inversion and rewrote the chain's
ISMIP7_MAP_OUT every iteration. The driver now refuses that combination before
it loads a mesh, so a check pointed at a live chain's MAP leaves it untouched.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DRIVER = REPO / "antarctica" / "scripts" / "inversion_icepack2.py"
REFUSAL = "ISMIP7_GRAD_CHECK=1 needs the TAO path"


def _run(tmp_path, **env):
    full = dict(os.environ, OMP_NUM_THREADS="1",
                ISMIP7_MESH=str(tmp_path / "no_such_mesh.msh"), **env)
    return subprocess.run([sys.executable, str(DRIVER)], cwd=REPO / "antarctica",
                          env=full, capture_output=True, text=True, timeout=300)


@pytest.mark.parametrize("precond", [None, "none", "mass"])
def test_grad_check_on_a_scipy_metric_is_refused_and_leaves_map_out(tmp_path, precond):
    map_out = tmp_path / "chain.h5"
    map_out.write_bytes(b"chain checkpoint")
    env = {"ISMIP7_GRAD_CHECK": "1", "ISMIP7_MAP_OUT": str(map_out),
           "ISMIP7_WARM_START": str(map_out)}
    if precond is not None:
        env["ISMIP7_GRAD_PRECOND"] = precond
    res = _run(tmp_path, **env)
    assert res.returncode != 0
    assert REFUSAL in res.stderr
    assert "Loading mesh" not in res.stdout
    assert map_out.read_bytes() == b"chain checkpoint"
    assert not (tmp_path / "chain.h5.done").exists()


def test_grad_check_on_the_tao_path_is_accepted(tmp_path):
    res = _run(tmp_path, ISMIP7_GRAD_CHECK="1", ISMIP7_GRAD_PRECOND="mass_consistent")
    # Past the knob checks: it fails on the missing mesh instead.
    assert REFUSAL not in res.stderr
    assert "Loading mesh" in res.stdout
