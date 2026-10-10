r"""Under ISMIP7_FLUIDITY_CONTROL=floating the inversion holds phi at zero on
the nodes where the grounded gate (1 - He) has shut it off
(``dual_friction.floating_control_nodes``), so the MAP carries no grounded
fluidity variation (issue #153)."""
import numpy as np
import firedrake as fd

from icepack2_tools.dual_friction import (
    GL_WIDTH, PHI_GROUNDED_TOL, floating_control_nodes)

RHO = 1024.0 / 917.0


def _strip(haf_left, haf_right, n=40):
    r"""A strip whose height above flotation runs linearly in x, per cell."""
    mesh = fd.RectangleMesh(n, 2, 40e3, 2e3)
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    x = fd.SpatialCoordinate(mesh)[0]
    H = fd.Function(Q0).interpolate(fd.Constant(500.0))
    haf = haf_left + (haf_right - haf_left) * x / 40e3
    b = fd.Function(Q0).interpolate(-(500.0 - haf) / RHO)
    xn = fd.Function(Q).interpolate(x).dat.data_ro.copy()
    return H, b, Q, xn, haf_left, haf_right


def test_grounded_nodes_are_held_and_the_band_stays_free():
    H, b, Q, xn, lo, hi = _strip(200.0, -200.0)
    free = floating_control_nodes(H, b, Q).dat.data_ro
    haf_n = lo + (hi - lo) * xn / 40e3
    # one cell is 1 km, 10 m of height above flotation: a node is held when
    # every cell it touches is past the tolerance (about 35 m at GL_WIDTH)
    cut = GL_WIDTH * np.arctanh(1.0 - 2.0 * PHI_GROUNDED_TOL)
    assert np.all(free[haf_n > cut + 10.0] == 0.0)
    assert np.all(free[haf_n < cut - 10.0] == 1.0)
    assert set(np.unique(free)) == {0.0, 1.0}


def test_all_floating_is_all_free_and_deep_grounding_is_all_held():
    H, b, Q, *_ = _strip(-50.0, -300.0)
    assert np.all(floating_control_nodes(H, b, Q).dat.data_ro == 1.0)
    H, b, Q, *_ = _strip(300.0, 100.0)
    assert np.all(floating_control_nodes(H, b, Q).dat.data_ro == 0.0)
