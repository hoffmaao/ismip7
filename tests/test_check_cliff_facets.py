r"""The cliff-facet census (``antarctica/scripts/check_cliff_facets.py``,
issue #166) on a slab with a known answer: 500 m of ice on a 100 m bed in
x < 10 km, and beside it four 2 km bands of ice-free ground, lower land (bed
0 m), a partial wall (300 m), rock above the ice surface (700 m) and shallow
sea (-50 m). Each band meets the ice along two 1 km facets."""
import importlib.util
import os
import sys

import numpy as np
import pytest

fd = pytest.importorskip("firedrake")

from icepack2.constants import gravity as g, ice_density as rho_I  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "antarctica", "scripts")
H_ICE, B_ICE, X_FRONT = 500.0, 100.0, 10e3
BANDS = (("land_low", 0.0), ("land_partial", 300.0), ("land_above", 700.0), ("ocean_deep", -50.0))


def _script():
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    spec = importlib.util.spec_from_file_location(
        "check_cliff_facets", os.path.join(SCRIPTS, "check_cliff_facets.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _state(path, speed):
    mesh = fd.RectangleMesh(16, 8, 16e3, 8e3)
    x, y = fd.SpatialCoordinate(mesh)
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    ice = fd.lt(x, X_FRONT)
    bed = fd.conditional(ice, B_ICE, fd.conditional(
        fd.lt(y, 2e3), BANDS[0][1], fd.conditional(
            fd.lt(y, 4e3), BANDS[1][1], fd.conditional(fd.lt(y, 6e3), BANDS[2][1], BANDS[3][1]))))
    b = fd.Function(Q0, name="bed").interpolate(bed)
    H = fd.Function(Q0, name="thickness").interpolate(fd.conditional(ice, H_ICE, 0.0))
    s = fd.Function(Q0, name="surface").interpolate(fd.max_value(b + H, 0.1 * H))
    u = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1), name="velocity").interpolate(
        fd.as_vector((speed, 0.0)))
    tau = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0), name="basal_stress").interpolate(
        fd.as_vector((-0.05, 0.0)))                                  # 50 kPa
    with fd.CheckpointFile(str(path), "w") as chk:
        chk.save_mesh(mesh)
        for f in (H, b, s, u, tau):
            chk.save_function(f)
    return str(path)


def _exposed(b_rock):
    s = B_ICE + H_ICE
    h_e = min(H_ICE, max(s - b_rock, 0.0))
    return 0.5 * g * rho_I * h_e ** 2


def test_the_census_classes_each_band_and_gives_version_2_its_push(tmp_path):
    mod = _script()
    mesh, f = mod._load(_state(tmp_path / "slab.h5", 100.0))
    fc = mod.facet_census(mesh, f, 1.0)
    classes = mod.classify(fc["B"], fc["s"] - fc["H"], fc["s"])
    full = _exposed(-1e9)
    for name, b_rock in BANDS:
        m = classes[name]
        assert m.sum() == 2, name
        assert fc["length"][m].sum() == pytest.approx(2e3)
        assert np.allclose(fc["push1"][m], full, rtol=1e-12)
        assert np.allclose(fc["push2"][m], _exposed(b_rock), rtol=1e-12, atol=1e-12 * full)
    assert sum(m.sum() for m in classes.values()) == fc["length"].size == 8
    assert fc["push2"][classes["land_partial"]][0] / full == pytest.approx(0.36)


def test_the_edge_cells_compare_the_change_with_their_drag(tmp_path):
    mod = _script()
    mesh, f = mod._load(_state(tmp_path / "slab.h5", 100.0))
    fc = mod.facet_census(mesh, f, 1.0)
    ec = mod.edge_cell_ratio(mesh, f, 1.0)
    change = np.abs(fc["push2"] - fc["push1"]) * fc["length"]
    assert ec["dpush"].sum() == pytest.approx(change.sum(), rel=1e-12)
    assert np.allclose(ec["drag"], 0.05 * ec["area"], rtol=1e-12)
    # the partial wall and the rock above change; lower land and sea do not
    assert ec["dpush"].size == 4


def test_compare_reports_the_speed_change_and_the_discharge(tmp_path):
    mod = _script()
    mesh, f = mod._load(_state(tmp_path / "ref.h5", 100.0))
    ref, q_ref, qb_ref = mod.cell_state(mesh, f, 1.0, None)
    mesh2, f2 = mod._load(_state(tmp_path / "other.h5", 120.0))
    oth, q_oth, _ = mod.cell_state(mesh2, f2, 1.0, None)
    # grounded ice leaves through the 8 km front at u h
    assert q_ref == pytest.approx(100.0 * H_ICE * 8e3 * 917.0 / 1e12, rel=1e-9)
    assert q_oth / q_ref == pytest.approx(1.2, rel=1e-9)
    assert qb_ref == {}
    fc = mod.facet_census(mesh, f, 1.0)
    _, changed = mod.census_lines(fc, None)
    lines = mod.compare_lines(ref, oth, fc, changed, q_ref, q_oth, {}, {})
    assert "+20.000 %" in lines[0]
