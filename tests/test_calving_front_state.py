r"""A calving law reads its fields off a front state; the forward builds it
live (``simulation.calving_front_state``, an ``icepack_tools.calving.FrontState``)
and hands the law's rate to the shared level set as the prescribed rate."""
import importlib.util
import os

import numpy as np
from firedrake import (Constant, Function, FunctionSpace, TensorFunctionSpace,
                       UnitSquareMesh, VectorFunctionSpace, FiniteElement,
                       SpatialCoordinate, as_vector, conditional)
from icepack_tools import calving

from icepack2_tools.levelset import LevelSet

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "antarctica", "scripts")


def _simulation():
    import sys
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    spec = importlib.util.spec_from_file_location(
        "simulation", os.path.join(SCRIPTS, "simulation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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
    return mesh, z, h, b, ls


def _xc(mesh):
    return Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro[:, 0]


def test_the_fields_are_the_forwards_own():
    sim = _simulation()
    mesh, z, h, b, ls = _state()
    A = Constant(20.0)
    m = sim.calving_front_state(z, h, b, ls, A, 3.0)
    assert isinstance(m, calving.FrontState)
    assert m.u is z.subfunctions[0] and m.M is z.subfunctions[1] and m.tau is z.subfunctions[2]
    assert m.h is h and m.b is b and m.nfront is ls.ghat
    assert m.A is A and m.n == 3.0
    haf = Function(m.Q0).interpolate(m.haf).dat.data_ro
    chi = Function(m.Q0).interpolate(m.chi_gr).dat.data_ro
    xc = _xc(mesh)
    # 300 m of ice on a 100 m deep bed is grounded; on 600 m it floats. The
    # flotation test is the forward's own, under its seawater of 1024
    assert np.all(chi[(xc < 0.25)] == 1.0)
    assert np.all(chi[(xc > 0.25) & (xc < 0.5)] == 0.0)
    assert np.allclose(haf[xc < 0.25], 300.0 - 1024.0 / 917.0 * 100.0)
    assert calving.densities(m) == (sim.rho_I, sim.rho_W)


def test_a_law_rate_is_a_cell_field_the_level_set_can_take():
    sim = _simulation()
    mesh, z, h, b, ls = _state()
    m = sim.calving_front_state(z, h, b, ls)
    # the minimum-thickness law at Hc = 300 m: the front holds, the rate
    # equals the arrival speed on every marine front cell
    law = calving.make("thickness", hc=300.0)
    c = Function(m.Q0).interpolate(law.rate(m, 0.0)).dat.data_ro
    front = ls.front_len.dat.data_ro > 0.0
    assert front.any()
    assert np.allclose(c[front], 1e-3, rtol=1e-6)
    ls.advance(0.05, z.subfunctions[0], rate=law.rate(m, 0.0))
    beyond, frac = ls.calving_masks()
    assert frac is not None and frac.max() > 0.0
    assert not beyond.any()


def test_a_law_reads_the_front_normal_of_the_current_extent_at_the_first_advance():
    sim = _simulation()
    mesh, z, h, b, ls = _state()
    m = sim.calving_front_state(z, h, b, ls)
    ls.advance(0.05, z.subfunctions[0], rate=calving.make("velocity").rate(m, 0.0))
    front = ls.front_len.dat.data_ro > 0.0
    assert front.any()
    # the front at x = 0.5 faces +x and the ice flows at 1e-3 in +x
    assert np.allclose(ls.c_cell.dat.data_ro[front], 1e-3, rtol=0.05)
