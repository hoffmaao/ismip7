r"""The thermal fluidity prior over water: the base of floating ice sits at the
ice-ocean interface's melting point and supplies nothing more. No frictional
heat (the shelf's basal shear stress is zero), no geothermal flux, and no water
content above the melting point.

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

from icepack.constants import heat_capacity as c_heat, ice_density as rho_I  # noqa: E402

from icepack2_tools.thermo_model import (                       # noqa: E402
    DEFAULTS,
    OCEAN_T_FREEZE_0,
    OCEAN_T_FREEZE_GRAD,
    cap_floating_enthalpy,
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


def _energy(floating, water=2000.0, C=0.0, q_geo=0.0, speed=300.0, acc=0.1):
    mesh, Q, Q0, h, bed, s, u = _column(floating, water)
    u.interpolate(Constant((speed, 0.0)))
    p = dict(DEFAULTS, q_geo=q_geo)
    return solve_energy(u, h, s, Function(Q).assign(C), Function(Q0).assign(10.0), bed,
                        Function(Q0).assign(acc), p, with_strain=False,
                        T_srf_field=None).dat.data_ro.copy()


def _heat_from_friction(floating, water=2000.0):
    r"""Enthalpy with a friction coefficient minus enthalpy without one."""
    return _energy(floating, water, C=0.05) - _energy(floating, water, C=0.0)


def test_friction_warms_a_grounded_column():
    assert np.all(_heat_from_friction(floating=False) > 0.0)


def test_friction_does_not_warm_ice_over_water():
    r"""Over deep water, and over the thin water column of a grounding zone,
    where the old gate still passed 40 % of the frictional heat."""
    for water in (2000.0, 10.0):
        assert np.all(_heat_from_friction(floating=True, water=water) == 0.0)


def test_geothermal_flux_warms_the_bed_and_not_ice_over_water():
    for floating, warmer in ((False, True), (True, False)):
        for water in ((2000.0, 10.0) if floating else (2000.0,)):
            dE = _energy(floating, water, q_geo=50.0) - _energy(floating, water, q_geo=0.0)
            assert np.all(dE > 0.0) if warmer else np.all(dE == 0.0)


def test_a_floating_column_settles_between_its_surface_and_the_ocean():
    r"""No sources, no flow, no accumulation: the surface and the ice-ocean
    interface pull equally, so the column sits midway between them, while a
    grounded column, with no basal exchange, sits at its surface value."""
    T_srf = DEFAULTS["T_srf"]
    T_f = OCEAN_T_FREEZE_0 - OCEAN_T_FREEZE_GRAD * RATIO * H_ICE
    for floating, want in ((True, 0.5 * (T_srf + T_f)), (False, T_srf)):
        E = _energy(floating, speed=0.0, acc=0.0)
        assert np.allclose(E / float(rho_I * c_heat), want, rtol=1e-10)


def test_floating_ice_holds_no_water():
    mesh, Q, Q0, h, bed, s, u = _column(True)
    E_melt = float(rho_I * c_heat) * 273.15
    E = Function(Q0).assign(1.02 * E_melt)          # temperate with water
    g = Function(Q0).assign(0.0)
    g.dat.data[: len(g.dat.data) // 2] = 1.0        # half the cells on the bed
    cap_floating_enthalpy(E, g)
    on_bed = g.dat.data_ro > 0.5
    assert np.allclose(E.dat.data_ro[~on_bed], E_melt)
    assert np.allclose(E.dat.data_ro[on_bed], 1.02 * E_melt)
