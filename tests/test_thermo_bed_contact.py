r"""The thermal fluidity prior heats with friction only where the ice base
rests on the bed: ice sliding over water does no frictional work on itself.

The earlier gate, ``1 - 0.5 (1 + tanh(gap / 50 m))``, was centred on a zero
gap, where every grounded column sits under the surface this model builds, so
grounded ice took half its frictional heat and floating ice near the grounding
line took a share. Serial, a 40 km x 20 km rectangle, no data files.
"""

# icepack2 -> irksome first, before any UFL form is assembled
import icepack2_tools.dual_friction  # noqa: F401,E402

import numpy as np                                              # noqa: E402

from firedrake import (                                         # noqa: E402
    Constant,
    Function,
    FunctionSpace,
    RectangleMesh,
    VectorFunctionSpace,
)

from icepack2_tools.thermo_model import (                       # noqa: E402
    DEFAULTS,
    grounded_frac,
    solve_energy,
)

RATIO = 917.0 / 1024.0
H_ICE = 800.0


def _column(floating, water=2000.0):
    r"""An 800 m column, grounded on a bed at -100 m, or floating over
    ``water`` metres of ocean."""
    mesh = RectangleMesh(20, 10, 40e3, 20e3)
    Q = FunctionSpace(mesh, "CG", 1)
    Q0 = FunctionSpace(mesh, "DG", 0)
    h = Function(Q).assign(H_ICE)
    bed = Function(Q).assign(-(RATIO * H_ICE + water) if floating else -100.0)
    s = Function(Q).assign((1 - RATIO) * H_ICE if floating else -100.0 + H_ICE)
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(Constant((300.0, 0.0)))
    return mesh, Q, Q0, h, bed, s, u


def test_the_gate_is_one_on_the_bed_and_zero_over_water():
    for floating, want in ((False, 1.0), (True, 0.0)):
        mesh, Q, Q0, h, bed, s, u = _column(floating)
        g = Function(Q0).interpolate(grounded_frac(h, s, bed)).dat.data_ro
        assert np.all(g == want)


def _heat_from_friction(floating, water=2000.0):
    r"""Enthalpy with a friction coefficient minus enthalpy without one."""
    mesh, Q, Q0, h, bed, s, u = _column(floating, water)
    A = Function(Q0).assign(10.0)
    acc = Function(Q0).assign(0.1)
    E = {}
    for c in (0.0, 0.05):
        C = Function(Q).assign(c)
        E[c] = solve_energy(u, h, s, C, A, bed, acc, dict(DEFAULTS),
                            with_strain=False, T_srf_field=None).dat.data_ro.copy()
    return E[0.05] - E[0.0]


def test_friction_warms_a_grounded_column():
    assert np.all(_heat_from_friction(floating=False) > 0.0)


def test_friction_does_not_warm_ice_over_water():
    r"""Over deep water, and over the thin water column of a grounding zone,
    where the old gate still passed 40 % of the frictional heat."""
    for water in (2000.0, 10.0):
        assert np.all(_heat_from_friction(floating=True, water=water) == 0.0)
