r"""The relative-decrease stopping rule, and the L-surface driver's pure parts."""
import importlib.util
import json
import os

import numpy as np
import pytest

from icepack2_tools.optimization import FunctionalDecreaseStop

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _lsurface():
    spec = importlib.util.spec_from_file_location(
        "lsurface", os.path.join(REPO, "antarctica", "scripts", "lsurface.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── FunctionalDecreaseStop ───────────────────────────────────────────────

def test_the_first_iterate_never_stops():
    assert FunctionalDecreaseStop(1e-4).update(0, 100.0) is False


def test_the_criterion_is_the_relative_decrease():
    stop = FunctionalDecreaseStop(ftol=1e-4, min_iter=0)
    stop.update(0, 200.0)
    stop.update(1, 100.0)
    assert stop.criterion == pytest.approx((200.0 - 100.0) / 200.0)


def test_small_functionals_are_measured_against_one_not_themselves():
    """max(|J_old|, |J_new|, 1): a functional already below 1 is judged on
    its absolute decrease, as in scipy."""
    stop = FunctionalDecreaseStop(ftol=1e-4, min_iter=0)
    stop.update(0, 0.02)
    stop.update(1, 0.01)
    assert stop.criterion == pytest.approx(0.01)


def test_it_stops_when_the_decrease_falls_below_ftol():
    stop = FunctionalDecreaseStop(ftol=1e-3, min_iter=0)
    stop.update(0, 1000.0)
    assert stop.update(1, 990.0) is False        # 1e-2 > ftol
    assert stop.update(2, 989.5) is True         # ~5e-4 <= ftol


def test_min_iter_holds_it_back():
    stop = FunctionalDecreaseStop(ftol=1e-3, min_iter=3)
    stop.update(0, 1000.0)
    assert stop.update(1, 1000.0) is False       # flat, but it < min_iter
    assert stop.update(2, 1000.0) is False
    assert stop.update(3, 1000.0) is True


def test_zero_disables_the_rule():
    stop = FunctionalDecreaseStop(ftol=0.0, min_iter=0)
    stop.update(0, 1.0)
    assert stop.update(1, 1.0) is False


# ── lsurface ─────────────────────────────────────────────────────────────

def test_the_submit_line_carries_both_weights_the_tolerance_and_distinct_outputs(tmp_path):
    ls = _lsurface()
    cmd = ls.submit_command("budd", 1e4, 1e6, str(tmp_path), "1e-4",
                            ["ISMIP7_LC=2000", "ISMIP7_MESH=/m.msh"], dry_run=True)
    assert cmd[:3] == [ls.SUBMIT, "inversion", "--name"]
    assert cmd[3] == "ls_budd_t4_p6"
    assert "--dry-run" in cmd
    kv = dict(a.split("=", 1) for a in cmd if "=" in a)
    assert kv["ISMIP7_FRICTION"] == "budd"
    assert kv["ISMIP7_GAMMA_THETA"] == "10000"
    assert kv["ISMIP7_GAMMA_PHI"] == "1e+06"
    assert kv["ISMIP7_FTOL"] == "1e-4"
    assert kv["ISMIP7_MAP_OUT"] == str(tmp_path / "gt10000_gp1e+06.h5")
    assert kv["ISMIP7_INVERSION_TIMING_JSON"] == str(tmp_path / "gt10000_gp1e+06.json")
    assert kv["ISMIP7_LC"] == "2000" and kv["ISMIP7_MESH"] == "/m.msh"


def test_a_relative_out_dir_reaches_the_job_as_an_absolute_path(tmp_path, monkeypatch):
    """submit.sh changes directory before sbatch, so a relative path would
    resolve against the repo root in the job and the caller's cwd in harvest."""
    ls = _lsurface()
    calls = []
    monkeypatch.setattr(ls.subprocess, "run", lambda cmd, check: calls.append(cmd))
    monkeypatch.chdir(tmp_path)
    ls.main(["submit", "--gammas", "1e3", "--out-dir", "rel/x"])
    (cmd,) = calls
    kv = dict(a.split("=", 1) for a in cmd if "=" in a)
    assert kv["ISMIP7_MAP_OUT"] == str(tmp_path / "rel" / "x" / "gt1000_gp1000.h5")
    assert kv["ISMIP7_INVERSION_TIMING_JSON"] == str(tmp_path / "rel" / "x" / "gt1000_gp1000.json")
    assert (tmp_path / "rel" / "x").is_dir()


def test_two_points_never_share_a_map_or_a_json(tmp_path):
    ls = _lsurface()
    a = ls.point_paths(str(tmp_path), 1e3, 1e5)
    b = ls.point_paths(str(tmp_path), 1e5, 1e3)
    assert len({*a, *b}) == 4


def test_the_corner_is_where_two_power_laws_meet():
    """A misfit that is flat for small gamma and then rises, and a prior
    energy that falls steeply and then flattens, meet at gamma = 1: that joint
    is the point of maximum curvature, whichever way the points are ordered."""
    ls = _lsurface()
    gammas = np.logspace(-3, 3, 13)
    misfit = 1.0 + gammas ** 2                 # flat, then rises steeply
    reg = 1.0 + 1.0 / gammas ** 2              # falls steeply, then flat
    k = ls.lcurve_corner(gammas, misfit, reg)
    assert abs(np.log10(gammas[k])) <= 0.5
    shuffled = np.random.default_rng(0).permutation(len(gammas))
    k2 = ls.lcurve_corner(gammas[shuffled], misfit[shuffled], reg[shuffled])
    assert np.isclose(gammas[shuffled][k2], gammas[k])


def test_a_concave_kink_is_not_a_corner():
    """The same curve with the axes swapped turns the other way; a signed
    curvature finds no convex corner there to prefer over the ends."""
    ls = _lsurface()
    gammas = np.logspace(-3, 3, 13)
    misfit = 1.0 + gammas ** 2
    reg = 1.0 + 1.0 / gammas ** 2
    k = ls.lcurve_corner(gammas, reg, misfit)
    assert abs(np.log10(gammas[k])) > 0.5


def test_coincident_points_cannot_win_the_corner():
    ls = _lsurface()
    assert ls.lcurve_corner([1, 10, 100], [2, 2, 2], [5, 5, 5]) == -1


def test_fewer_than_three_points_have_no_corner():
    ls = _lsurface()
    assert ls.lcurve_corner([1, 10], [2, 1], [1, 2]) == -1


def test_harvest_reads_the_inversions_own_json_records(tmp_path):
    """load_points consumes what inversion_icepack2.py writes: the knobs block
    for the weights and the last evaluation for the terms."""
    ls = _lsurface()
    for gt, gp, mis in ((1e3, 1e5, 5.0), (1e5, 1e5, 9.0)):
        _, fn = ls.point_paths(str(tmp_path), gt, gp)
        with open(fn, "w") as f:
            json.dump({
                "phase": "finished", "nit": 7, "message": "CONVERGED",
                "knobs": {"gamma_theta": gt, "gamma_phi": gp},
                "evaluations": [
                    {"eval": 0, "misfit": 50.0, "reg_theta": 0.0, "reg_phi": 0.0},
                    {"eval": 7, "misfit": mis, "reg_theta": 1.0, "reg_phi": 2.0},
                ],
            }, f)
    rows = ls.load_points(str(tmp_path))
    assert [(r["gamma_theta"], r["misfit"], r["nit"], r["finished"]) for r in rows] == [
        (1e3, 5.0, 7, True), (1e5, 9.0, 7, True)]
    axes = [(axis, fixed, k) for axis, fixed, *_, k in ls.slices(rows)]
    # one gamma_phi row with two gamma_theta points: a slice, but no corner yet
    assert ("theta", 1e5, -1) in axes


def _write_record(ls, out_dir, gt, gp, mis, rt, rp, phase):
    _, fn = ls.point_paths(str(out_dir), gt, gp)
    with open(fn, "w") as f:
        json.dump({
            "phase": phase, "nit": 7, "message": "",
            "knobs": {"gamma_theta": gt, "gamma_phi": gp},
            "evaluations": [{"eval": 7, "misfit": mis, "reg_theta": rt, "reg_phi": rp}],
        }, f)


def test_harvest_lists_unfinished_points_but_fits_only_finished_ones(tmp_path, capsys):
    ls = _lsurface()
    _write_record(ls, tmp_path, 1e3, 1e5, 5.0, 3.0, 1.0, "finished")
    _write_record(ls, tmp_path, 1e4, 1e5, 9.0, 2.0, 1.0, "running")
    rows = ls.load_points(str(tmp_path))
    assert [(r["gamma_theta"], r["phase"], r["finished"]) for r in rows] == [
        (1e3, "finished", True), (1e4, "running", False)]
    ls.main(["harvest", str(tmp_path)])
    out = capsys.readouterr().out
    assert "skipped 1 unfinished point(s)" in out
    assert "running" in out
    with open(tmp_path / "lsurface.csv") as f:
        assert [line.rsplit(",", 1)[1].strip() for line in f][1:] == ["finished", "running"]


def test_harvest_with_no_finished_point_says_so(tmp_path):
    ls = _lsurface()
    _write_record(ls, tmp_path, 1e3, 1e5, 5.0, 3.0, 1.0, "running")
    with pytest.raises(SystemExit, match="no finished points"):
        ls.main(["harvest", str(tmp_path)])
