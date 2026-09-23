r"""A diagnostic solve that converges onto a speed no ice can have must be
treated as a failed solve, so the step reaches the rescue ladder and the
subcycles instead of being carried forward.

Measured motivation: a 1 km control in the Lambert/Amery grounding trough
reported SNES reason=2 with a function norm below 1 and a peak speed of
2.5e6 m/yr, and the step loop accepted it, because it asked only whether the
solver had converged.

Serial, one small rectangle mesh, no data files.
"""
import pytest

from icepack2_tools.runconfig import U_LIM_DEFAULT, max_speed_bound

fd = pytest.importorskip("firedrake")

from icepack2_tools.speed_bound import (                        # noqa: E402
    limiter_speed_bound,
    speed_bound_violation,
)


def test_the_default_bound_is_above_any_real_ice_and_below_a_runaway(monkeypatch):
    monkeypatch.delenv("ISMIP7_MAX_SPEED", raising=False)
    bound = max_speed_bound()
    assert bound == pytest.approx(20000.0)
    # the fastest ice measured anywhere is about 17 km/yr
    assert bound > 17000.0
    # the trough ignited at 20435 m/yr, one step before the solver failed
    assert bound < 20435.0


def test_the_bound_can_be_disabled_but_not_inverted(monkeypatch):
    monkeypatch.setenv("ISMIP7_MAX_SPEED", "0")
    assert max_speed_bound() == 0.0
    monkeypatch.setenv("ISMIP7_MAX_SPEED", "-1")
    with pytest.raises(ValueError, match="positive"):
        max_speed_bound()
    monkeypatch.setenv("ISMIP7_MAX_SPEED", "5e4")
    assert max_speed_bound() == pytest.approx(5e4)


def test_a_limiter_rung_is_not_vetoed_where_the_limiter_pins_a_node(monkeypatch):
    r"""A pinned node settles at 2.01e4 to 2.1e4 m/yr with the default u_lim,
    above the default bound; the limiter rung must accept it and still reject
    the runaway."""
    monkeypatch.delenv("ISMIP7_MAX_SPEED", raising=False)
    bound = limiter_speed_bound(max_speed_bound(), float(U_LIM_DEFAULT))
    assert bound > 2.1e4
    assert bound < 2.48e6
    assert limiter_speed_bound(1e5, 2e4) == pytest.approx(1e5)
    assert limiter_speed_bound(0.0, 2e4) == 0.0


def _field(speed_x):
    mesh = fd.RectangleMesh(8, 4, 8e4, 4e4)
    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    x, y = fd.SpatialCoordinate(mesh)
    u = fd.Function(V).interpolate(fd.as_vector((speed_x(x), 0.0 * y)))
    speed = fd.Function(Q)
    xy = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        fd.SpatialCoordinate(mesh))
    return u, speed, xy


def test_a_runaway_is_rejected_with_its_location():
    # 5 km/yr everywhere except a trough node at x = 6e4 running away
    u, speed, xy = _field(
        lambda x: fd.conditional(abs(x - 6e4) < 1.0, 2.5e6, 5e3))
    violation = speed_bound_violation(u, speed, xy, 2e4)
    assert violation is not None
    u_max, (x_max, _y_max) = violation
    assert u_max == pytest.approx(2.5e6)
    assert x_max == pytest.approx(6e4)


def test_a_healthy_solve_is_accepted():
    u, speed, xy = _field(lambda x: 5.5e3 * x / 8e4)
    assert speed_bound_violation(u, speed, xy, 2e4) is None


def test_the_limiter_bound_accepts_a_pinned_node_the_ordinary_bound_rejects():
    u, speed, xy = _field(
        lambda x: fd.conditional(abs(x - 6e4) < 1.0, 2.05e4, 5e3))
    assert speed_bound_violation(u, speed, xy, 2e4) is not None
    assert speed_bound_violation(
        u, speed, xy, limiter_speed_bound(2e4, 2e4)) is None


def test_a_disabled_bound_accepts_anything_without_touching_the_field():
    u, speed, xy = _field(lambda x: 2.5e6 + 0.0 * x)
    speed.assign(-1.0)
    assert speed_bound_violation(u, speed, xy, 0.0) is None
    assert speed.dat.data_ro.max() == -1.0
