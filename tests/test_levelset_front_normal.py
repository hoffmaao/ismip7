r"""A ``prescribed`` rate that reads the level set's front normal ``ghat``,
as the horizontal-force-balance law does, sees the normal of the current
extent, including at the first advance."""
import numpy as np
from firedrake import (Constant, Function, FunctionSpace, UnitSquareMesh,
                       VectorFunctionSpace, SpatialCoordinate, as_vector, dot,
                       conditional, max_value)

from icepack2_tools.levelset import LevelSet


def test_a_rate_reads_the_front_normal_of_the_current_extent_at_the_first_advance():
    mesh = UnitSquareMesh(8, 8)
    x, _ = SpatialCoordinate(mesh)
    # a speed small against the unit domain, so one 0.05 yr step retreats the
    # front by a fraction of a cell rather than past the whole ice body
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((Constant(1e-3), Constant(0.0))))
    Q0 = FunctionSpace(mesh, "DG", 0)
    # ice on the left half, 300 m thick; ocean on the right
    h = Function(Q0).interpolate(conditional(x < 0.5, Constant(300.0), Constant(0.0)))
    b = Function(Q0).interpolate(Constant(-600.0))
    ls = LevelSet(mesh, h, law="prescribed", h_min=1.0, anchor="extent")
    rate = max_value(dot(u, ls.ghat), Constant(0.0))
    ls.advance(0.05, u, h, b, rate=rate)
    front = ls.front_len.dat.data_ro > 0.0
    assert front.any()
    # the front at x = 0.5 faces +x and the ice flows at 1e-3 in +x
    assert np.allclose(ls.c_cell.dat.data_ro[front], 1e-3, rtol=0.05)
