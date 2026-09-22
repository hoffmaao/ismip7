r"""A hoffmaao/calving law reads seven fields from its model; the forward
supplies them live (simulation.LiveCalvingState) and hands the law's rate
to the shared level set as the prescribed rate."""
import importlib.util
import os

import numpy as np
from firedrake import (Constant, Function, FunctionSpace, TensorFunctionSpace,
                       UnitSquareMesh, VectorFunctionSpace, FiniteElement,
                       SpatialCoordinate, as_vector, dot, conditional, gt)

from icepack2_tools.levelset import LevelSet

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "antarctica", "scripts")


def _live_state_class():
    import sys
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    spec = importlib.util.spec_from_file_location(
        "simulation", os.path.join(SCRIPTS, "simulation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.LiveCalvingState


def _state():
    mesh = UnitSquareMesh(8, 8)
    x, y = SpatialCoordinate(mesh)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Sigma = TensorFunctionSpace(mesh, dg0, symmetry=True)
    T = VectorFunctionSpace(mesh, dg0)
    Z = V * Sigma * T
    z = Function(Z)
    # a speed small against the unit domain, so one 0.05 yr step retreats the
    # front by a fraction of a cell rather than past the whole ice body
    z.sub(0).interpolate(as_vector((Constant(1e-3), Constant(0.0))))
    Q0 = FunctionSpace(mesh, "DG", 0)
    # ice on the left half, 300 m thick; ocean on the right
    h = Function(Q0).interpolate(conditional(x < 0.5, Constant(300.0), Constant(0.0)))
    b = Function(Q0).interpolate(conditional(x < 0.25, Constant(-100.0), Constant(-600.0)))
    ls = LevelSet(mesh, h, law="prescribed", h_min=1.0, anchor="extent")
    ls._update_unit_gradient()
    return mesh, z, h, b, ls


def test_the_seven_fields_are_the_forwards_own():
    Live = _live_state_class()
    mesh, z, h, b, ls = _state()
    m = Live(z, h, b, ls)
    assert m.u is z.subfunctions[0] and m.M is z.subfunctions[1] and m.tau is z.subfunctions[2]
    assert m.h is h and m.nfront is ls.ghat
    haf = Function(m.Q0).interpolate(m.haf).dat.data_ro
    chi = Function(m.Q0).interpolate(m.chi_gr).dat.data_ro
    xc = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro[:, 0]
    # 300 m of ice on a 100 m deep bed is grounded; on 600 m it floats
    assert np.all(chi[(xc < 0.25)] == 1.0)
    assert np.all(chi[(xc > 0.25) & (xc < 0.5)] == 0.0)
    assert np.allclose(haf[xc < 0.25], 300.0 - 1028.0 / 917.0 * 100.0)


def test_a_law_rate_is_a_cell_field_the_level_set_can_take():
    Live = _live_state_class()
    mesh, z, h, b, ls = _state()
    m = Live(z, h, b, ls)

    class SpeedLaw:
        def rate(self, model, t):
            return dot(model.u, model.u) ** 0.5 * (Constant(1.0) - model.chi_gr)

        def describe(self):
            return "speed on floating cells"

    law = SpeedLaw()
    c = Function(m.Q0).interpolate(law.rate(m, 0.0)).dat.data_ro
    xc = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro[:, 0]
    assert np.allclose(c[(xc > 0.25) & (xc < 0.5)], 1e-3)
    assert np.all(c[xc < 0.25] == 0.0)
    ls.advance(0.05, z.subfunctions[0], h, b, rate=law.rate(m, 0.0))
    beyond, frac = ls.calving_masks()
    assert frac is not None and frac.max() > 0.0
    assert not beyond.any()
