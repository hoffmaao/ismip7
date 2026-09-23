r"""A calving result has to carry its threshold. The von Mises rate is
inversely proportional to sigma_max, so a run recorded as "vonmises" and
nothing else cannot be reproduced, and the default has changed before."""
import os
import re
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "antarctica", "scripts"))

import firedrake as fd                                           # noqa: E402

from icepack2_tools.runconfig import calving_sigma_max           # noqa: E402
from simulation import calving_front_owner, save_model_state     # noqa: E402

_SIGMA = re.compile(r"sigma_max grounded (\S+) MPa floating (\S+) MPa")


def test_von_mises_owner_round_trips_the_exact_thresholds(monkeypatch):
    r"""0.1234567 needs seven significant digits; `%g` keeps six and would
    record 0.123457, a different threshold from the one the run used."""
    monkeypatch.setenv("ISMIP7_CALVING_SIGMA_MAX_GROUNDED", "1.0000001")
    monkeypatch.setenv("ISMIP7_CALVING_SIGMA_MAX_FLOATING", "0.1234567")
    owner = calving_front_owner("vonmises", False, calving_sigma_max())
    grounded, floating = _SIGMA.search(owner).groups()
    assert float(grounded) == 1.0000001
    assert float(floating) == 0.1234567
    assert owner.startswith("level-set vonmises law (ISMIP7_CALVING=vonmises")


def test_owner_without_a_threshold_law():
    assert calving_front_owner("fixed", True) == (
        "level-set fixed law (ISMIP7_CALVING=fixed); "
        "ISMIP7_FIXED_FRONT is set but ignored for removal")
    assert calving_front_owner("none", True) == (
        "legacy fixed-front mask (ISMIP7_FIXED_FRONT)")
    assert calving_front_owner("none", False) == "none (no calving sink)"


def _minimal_ctx():
    mesh = fd.UnitSquareMesh(2, 2)
    Q = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    Z = V * fd.TensorFunctionSpace(mesh, "DG", 0) * fd.VectorFunctionSpace(
        mesh, "DG", 0)
    ctx = {"mesh": mesh, "z": fd.Function(Z), "geom_dg": True}
    for name in ("h", "theta", "phi", "b", "s", "phi_eff"):
        ctx[name] = fd.Function(Q)
    ctx["u_obs"] = fd.Function(V)
    return ctx


def test_the_checkpoint_carries_the_front_owner(tmp_path, monkeypatch):
    r"""The log is the first thing lost, so the data file must carry it too."""
    monkeypatch.setenv("ISMIP7_CALVING_SIGMA_MAX_FLOATING", "0.1234567")
    ctx = _minimal_ctx()
    ctx["calving_front_owner"] = calving_front_owner(
        "vonmises", False, calving_sigma_max())
    path = str(tmp_path / "state.h5")
    save_model_state(ctx, path, 0.0)
    with fd.CheckpointFile(path, "r") as chk:
        recorded = chk.get_attr("/", "calving_front_owner")
    assert recorded == ctx["calving_front_owner"]
    assert float(_SIGMA.search(recorded).group(2)) == 0.1234567
