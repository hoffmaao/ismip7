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

The term is written on ``div(u - u_ref)``, the surface change linearized
about a reference velocity ``u_ref``. Linearized about any ``u_ref``, the
surface a step ahead is

    s + theta tau gamma (T(u_ref) - h div(u - u_ref))

with ``T(u_ref)`` the thickness tendency the transport applies with
``u_ref``. The reference only sets the right-hand side: the Jacobian, and
with it the stability above, is the same for every ``u_ref``.
``ISMIP7_FSSA_REFERENCE`` picks it (:func:`resolve_reference`):

``start``
    ``u_ref`` is the velocity the run's balanced apparent mass balance was
    built from: there ``a + a_ref = div(u_ref h)`` cell by cell, so
    ``T(u_ref) = 0`` at t = 0, the dropped source and the divergence cancel
    exactly at the reference state, and the stabilized forward reproduces
    its MAP velocity at t = 0. Later the term omits ``T(u_ref)`` at the
    evolved thickness and forcing, an O(theta dt) load that stays at a new
    steady state. Without the apparent mass balance the starting state is
    out of balance, ``T(u_ref)`` is the full initial imbalance, and the term
    pulls ``div(u)`` toward the MAP's own divergence for the whole run.
``step``
    ``u_ref`` is the velocity the last advance used, and ``T`` that
    advance's realized tendency, ``(h_new - h_old) / dt``, which the
    implicit upwind transport makes ``source - div_FV(h_new u_ref)`` to its
    solver tolerance, with the source as the positivity limiter left it.
    The load vanishes at every steady state of the transport, whatever the
    forcing, so the long-run state does not move with dt or theta. Under a
    balanced apparent mass balance the first advance changes nothing, so
    ``T = 0`` and ``u_ref`` is the MAP velocity at the first solve, the
    same residual as ``start``.
``auto``
    ``start`` when the run carries an apparent mass balance, ``step`` when
    it does not.

Measured at 32 km without an apparent mass balance and without forcing, ten
years from the MAP's mixed state, against the unstabilized forward at dt
0.0125 (records ``test-32km-fssa-*``):

    ============  =====================  ======================  ========
    dt (yr)       thickness RMS (m)      worst cell (m)          mass (Gt)
    ============  =====================  ======================  ========
    start 0.05    0.47                   177                     -55
    step 0.05     0.06                   4.5                     +16
    start 0.1     0.91                   288                     -110
    step 0.1      0.11                   8.1                     +26
    start 0.2     1.84                   393                     -211
    step 0.2      0.20                   14.6                    +43
    ============  =====================  ======================  ========

The worst ``start`` cell is on the Antarctic Peninsula, at 788 m against the
reference's 395 m after ten years at dt 0.2; ``step`` leaves it at 399 m.
The unstabilized forward, stable at this resolution, misses by 0.07 and
0.15 m RMS at dt 0.1 and 0.2. With the apparent mass balance, ``auto``
reproduces ``start`` to the run-to-run noise of the same code (9e-13 m).
"""
from firedrake import Constant, TestFunction, div, dx, split

from icepack_tools.constants import gravity, ice_density, water_density
from icepack_tools.grounding import grounded_mask

REFERENCES = ("auto", "start", "step")
# Begins the line every forward prints with the weight it steps with
# (fssa_banner), which core_report lifts into the run's report.
FSSA_MARKER = "Free-surface stabilization:"


def fssa_banner(theta):
    r"""The startup line stating the stabilization weight ``theta`` a forward
    steps with, the value :func:`solverconfig.forward_fssa_theta` resolved."""
    if theta > 0:
        return (f"{FSSA_MARKER} theta {theta:g} on the lagged "
                "thickness-velocity coupling (ISMIP7_FSSA_THETA)")
    return f"{FSSA_MARKER} off, theta 0 (ISMIP7_FSSA_THETA)"


def resolve_reference(requested, apparent_mb):
    r"""The reference ``start`` or ``step`` an ``ISMIP7_FSSA_REFERENCE``
    value names, ``auto`` resolved by whether the run carries an apparent
    mass balance (``apparent_mb``)."""
    if requested not in REFERENCES:
        raise ValueError(
            f"ISMIP7_FSSA_REFERENCE must be one of {', '.join(REFERENCES)}, "
            f"not {requested!r}")
    if requested == "auto":
        return "start" if apparent_mb else "step"
    return requested


def restart_reference_error(metadata, resolved, source):
    r"""The message refusing a restart from ``source`` whose checkpoint was
    stepped from another reference than ``resolved``, or None.

    ``metadata`` holds the checkpoint's attributes. A checkpoint stepped
    under the stabilization records ``fssa_reference``; one that records only
    ``fssa_tau`` was written before the choice existed, under ``start``; one
    with neither was stepped without the stabilization, and one whose
    ``fssa_tau`` is 0 was never stepped (a prepared timing cache), and any
    reference may start from either."""
    tau = metadata.get("fssa_tau")
    if tau is not None and float(tau) == 0.0:
        return None
    was = metadata.get("fssa_reference")
    if was is None and tau is not None:
        was = "start"
    if was is None or str(was) == resolved:
        return None
    return (f"restart checkpoint {source} was stepped with the free-surface "
            f"stabilization measured from `{was}`, and ISMIP7_FSSA_REFERENCE "
            f"resolves to `{resolved}` here; set ISMIP7_FSSA_REFERENCE={was} "
            f"to continue the run")


def fssa_term(z, u_ref, tau, H, b, *, tendency=None, gl_width=10.0,
              rho_I=ice_density, rho_W=water_density, g=gravity):
    r"""The stabilization form, or 0 when ``tau`` or ``u_ref`` is None.

    ``tau`` is a Constant holding ``theta * dt`` (yr), so a forward can set it
    per substep and zero it for a steady solve; ``u_ref`` a Function on the
    velocity space. ``tendency`` (m/yr, a cell field) is ``T(u_ref)``, the
    thickness tendency at the reference velocity; None leaves it out, which
    is the ``start`` reference. It adds a load and leaves the Jacobian as it
    is."""
    if tau is None or u_ref is None:
        return 0
    u = split(z)[0]
    v = split(TestFunction(z.function_space()))[0]
    He = grounded_mask(H, b, gl_width=gl_width, rho_I=rho_I, rho_W=rho_W)
    gamma = He + (Constant(1.0) - He) * Constant(1.0 - rho_I / rho_W)
    form = -tau * Constant(rho_I * g) * gamma * H ** 2 * div(u - u_ref) * div(v) * dx
    if tendency is not None:
        form += tau * Constant(rho_I * g) * gamma * H * tendency * div(v) * dx
    return form
