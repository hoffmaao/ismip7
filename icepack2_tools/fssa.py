r"""Free-surface stabilization of the lagged thickness-velocity coupling.

The forward advances the thickness with the velocity solved at the previous
geometry. That coupling is explicit in the velocity's response to the surface,
and the response is a diffusion: a cell thicker than its neighbours sheds ice
to them at a rate proportional to the excess, and a cell thinner than them
draws ice in. Stepped explicitly past 2 over that rate the excess flips sign
every step and grows, which is the period-2 flip-flop of the 2 km forwards.
The rate is large where the bed carries little drag and the ice is thick, and
it grows with resolution, so a fixed step that is safe in one place and year
is not safe in another.

The free-surface stabilization algorithm (Kaus, Mühlhaus and Schmid 2010 for
geodynamics; Löfgren, Ahlkrona and Helanow 2022 for ice sheets; Tominec and
Ahlkrona 2024 for a vertically integrated model) keeps the lagged ordering and
instead makes the momentum problem aware of the surface it is about to
produce: the gravity load is evaluated at the surface a step ahead,
``s + theta tau ds/dt`` with ``ds/dt = gamma (a - div(u h))`` and ``gamma``
the surface-to-thickness ratio (1 grounded, ``1 - rho_i/rho_w`` afloat). In
the depth-integrated momentum balance the extra load is
``rho g h grad(theta tau gamma (a - div(u h)))``, and integrating it by parts
on a cell-wise geometry gives the symmetric, positive form

    theta tau rho g gamma h^2 div(u) div(v) dx

on the velocity block, with the source part dropped. The residual here is
written as ``div(h M) + tau_b - rho g h grad(s)`` tested with ``v`` (the
membrane term enters as ``-h M : eps(v)``), so the stabilization, a bulk
viscosity on the divergence, enters with the same minus sign. With it the diffusive
response per step becomes ``tau lambda / (1 + theta tau lambda)`` instead of
``tau lambda``, so for ``theta = 1`` the lagged step is stable for every tau.

The term is written on ``div(u - u_ref)`` with ``u_ref`` the velocity the
run's balanced apparent mass balance was built from: there ``a + a_ref =
div(u_ref h)`` cell by cell, so the dropped source and the divergence cancel
exactly at the reference state and the stabilized forward reproduces its MAP
velocity at t = 0 as before. Away from it the term is the implicit surface
change the step produces, which is what the method means.
"""
from firedrake import Constant, TestFunction, div, dx, split

from icepack_tools.constants import gravity, ice_density, water_density
from icepack_tools.grounding import grounded_mask


def fssa_term(z, u_ref, tau, H, b, *, gl_width=10.0, rho_I=ice_density,
              rho_W=water_density, g=gravity):
    r"""The stabilization form, or 0 when ``tau`` or ``u_ref`` is None.

    ``tau`` is a Constant holding ``theta * dt`` (yr), so a forward can set it
    per substep and zero it for a steady solve; ``u_ref`` a Function on the
    velocity space."""
    if tau is None or u_ref is None:
        return 0
    u = split(z)[0]
    v = split(TestFunction(z.function_space()))[0]
    He = grounded_mask(H, b, gl_width=gl_width, rho_I=rho_I, rho_W=rho_W)
    gamma = He + (Constant(1.0) - He) * Constant(1.0 - rho_I / rho_W)
    return -tau * Constant(rho_I * g) * gamma * H ** 2 * div(u - u_ref) * div(v) * dx
