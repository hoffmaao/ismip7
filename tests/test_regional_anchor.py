r"""The friction anchor's driving stress: local, or the grounded ice's
regional average (``ISMIP7_ANCHOR_LENGTH``).

The local anchor ``tau_d / |u_obs|^(1/3)`` vanishes with the surface slope, so
the prior gives an ice divide no friction. Averaged over the grounded ice
within about the anchor length it stays finite across the crest, it leaves a
uniform slope as it is, and floating ice beyond a grounding line does not
dilute it.

Serial, a 120 km x 40 km rectangle with 2 km cells, no data files.
"""

# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
from icepack2_tools.dual_friction import (                      # noqa: E402
    rebase_log_friction,
    regional_driving_stress,
    weertman_anchor,
)
from icepack2.constants import ice_density as rho_I, gravity as g  # noqa: E402

import numpy as np                                              # noqa: E402
import pytest                                                   # noqa: E402

from firedrake import (                                         # noqa: E402
    Constant,
    Function,
    FunctionSpace,
    RectangleMesh,
    SpatialCoordinate,
    VectorFunctionSpace,
    grad,
    inner,
    max_value,
    sqrt,
)

L_X, L_Y = 120e3, 40e3
X0 = 60e3                       # the crest of the dome, and the grounding line in the shelf case


def _fields(surface, thickness, bed, speed=0.0):
    mesh = RectangleMesh(60, 20, L_X, L_Y)
    Q = FunctionSpace(mesh, "CG", 1)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)
    s = Function(Q).interpolate(surface(x))
    H = Function(Q).interpolate(thickness(x))
    b = Function(Q).interpolate(bed(x))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        Constant((speed, 0.0)))
    xc = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro[:, 0]
    return mesh, Q, Q0, s, H, b, u, xc


def _dome():
    # a smooth divide: 3000 m at the crest, sloping 0.001 at 60 km out
    return _fields(lambda x: 3000.0 - 30.0 * ((x - X0) / 60e3) ** 2,
                   lambda x: 2000.0 + 0.0 * x, lambda x: 1000.0 - 30.0 * ((x - X0) / 60e3) ** 2)


def test_zero_length_is_the_local_balance():
    mesh, Q, Q0, s, H, b, u, xc = _dome()
    u.interpolate(Constant((7.0, 0.0)))
    got = weertman_anchor(H, s, u, 3.0, Q0).dat.data_ro
    want = Function(Q0).interpolate(
        rho_I * g * H * sqrt(inner(grad(s), grad(s)) + Constant(1e-12))
        / max_value(sqrt(inner(u, u)), Constant(1.0)) ** (1.0 / 3.0)).dat.data_ro
    assert np.array_equal(got, want)


def test_a_divide_keeps_the_regional_friction():
    mesh, Q, Q0, s, H, b, u, xc = _dome()
    local = weertman_anchor(H, s, u, 3.0, Q0).dat.data_ro
    regional = weertman_anchor(H, s, u, 3.0, Q0, length=20e3, b=b).dat.data_ro
    crest = np.abs(xc - X0) < 1.5e3
    flank = np.abs(np.abs(xc - X0) - 40e3) < 1.5e3
    # the local anchor collapses at the crest; the regional one keeps a
    # sizeable share of the flank value
    assert local[crest].max() < 0.05 * local[flank].mean()
    assert regional[crest].min() > 0.25 * regional[flank].mean()


def test_a_uniform_slope_is_left_as_it_is():
    mesh, Q, Q0, s, H, b, u, xc = _fields(
        lambda x: 3000.0 - 0.002 * x, lambda x: 2000.0 + 0.0 * x, lambda x: 1000.0 - 0.002 * x)
    local = weertman_anchor(H, s, u, 3.0, Q0).dat.data_ro
    regional = weertman_anchor(H, s, u, 3.0, Q0, length=20e3, b=b).dat.data_ro
    assert np.allclose(regional, local, rtol=1e-6)


def test_floating_ice_does_not_dilute_the_grounded_average():
    r"""Grounded ice thinning onto a flat bed 500 m below sea level, reaching
    flotation at x = 60 km, and a flat shelf beyond with a continuous surface.
    The shelf has no driving stress; averaged in, it would pull the grounded
    cells near the grounding line to about half their value."""
    ratio = 917.0 / 1024.0
    h_float = 500.0 / ratio + 2.0                  # just above flotation on a -500 m bed
    mesh, Q, Q0, s, H, b, u, xc = _fields(
        lambda x: 0.0 * x, lambda x: 0.0 * x, lambda x: 0.0 * x)
    x, _ = SpatialCoordinate(mesh)
    from firedrake import conditional, lt
    grounded_x = lt(x, X0)
    H.interpolate(conditional(grounded_x, h_float + 0.01 * (X0 - x), h_float))
    b.interpolate(conditional(grounded_x, -500.0, -2000.0))
    s.interpolate(conditional(grounded_x, -500.0 + h_float + 0.01 * (X0 - x), (1.0 - ratio) * h_float))
    tau_c = Function(Q0).interpolate(regional_driving_stress(H, s, b, Q0, 20e3)).dat.data_ro
    local = Function(Q0).interpolate(
        rho_I * g * H * sqrt(inner(grad(s), grad(s)) + Constant(1e-12))).dat.data_ro
    near_gl = (xc < X0 - 2e3) & (xc > X0 - 8e3)
    assert np.all(tau_c[near_gl] > 0.9 * local[near_gl])


def test_the_anchor_length_is_checked():
    mesh, Q, Q0, s, H, b, u, xc = _dome()
    with pytest.raises(ValueError, match="needs the bed"):
        weertman_anchor(H, s, u, 3.0, Q0, length=20e3)
    with pytest.raises(ValueError, match=">= 0"):
        weertman_anchor(H, s, u, 3.0, Q0, length=-1.0, b=b)


def test_rebasing_theta_keeps_the_friction_on_grounded_ice():
    r"""ISMIP7_WARM_START_THETA=physical: a warm start's theta moved onto a new
    anchor gives back the same friction C exp(theta) on grounded ice, and
    leaves theta alone where the ice floats."""
    from firedrake import conditional, lt, sin
    mesh, Q, Q0, s, H, b, u, xc = _dome()
    x, _ = SpatialCoordinate(mesh)
    # grounded for x < 60 km, floating beyond (bed far below flotation)
    b.interpolate(conditional(lt(x, X0), 1000.0, -5000.0))
    theta = Function(Q).interpolate(0.3 * sin(x / 7e3))
    C_to = Function(Q0).assign(0.02)
    C_from = Function(Q0).assign(0.05)                 # the old anchor, 2.5x the new
    new = rebase_log_friction(theta, C_from, C_to, H, b)
    xn = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro[:, 0]
    inside = xn < X0 - 3e3              # nodes whose every cell is grounded
    shelf = xn > X0 + 3e3               # nodes whose every cell floats
    assert np.allclose(new.dat.data_ro[inside] - theta.dat.data_ro[inside], np.log(2.5))
    assert np.array_equal(new.dat.data_ro[shelf], theta.dat.data_ro[shelf])
