"""
Smooth Heaviside grounding zone treatment (Ua/Gudmundsson approach).

Provides functions for:
  - Smooth Heaviside of height-above-flotation
  - Blended surface elevation (grounded ↔ floating)
  - Effective pressure sliding coefficient (phi_eff)
  - Grounding line indicator for flux computation
"""

import firedrake as fd
from firedrake import Constant, Function, max_value, conditional


def height_above_flotation(H, b, rho_I=917.0, rho_W=1024.0):
    """Compute height above flotation: positive = grounded, negative = floating."""
    h_f = Constant(rho_W / rho_I) * max_value(-b, Constant(0.0))
    return H - h_f


def smooth_heaviside(haf, kH=1.0):
    """Smooth Heaviside of height-above-flotation.

    He = 0.5 + 0.5 * tanh(kH * haf)

    Parameters
    ----------
    haf : UFL expression
        Height above flotation (positive = grounded).
    kH : float
        Sharpness parameter (m⁻¹). 1/kH is the vertical transition width.

    Returns
    -------
    UFL expression for He ∈ [0, 1].
    """
    return Constant(0.5) + Constant(0.5) * fd.tanh(Constant(kH) * haf)


def blended_surface(H, b, He, rho_ratio=917.0 / 1024.0):
    """Surface elevation blended between grounded and floating.

    s = He * (b + H) + (1 - He) * (1 - rho_ratio) * H

    Parameters
    ----------
    H : Function or UFL expression
        Ice thickness.
    b : Function
        Bed elevation.
    He : UFL expression
        Smooth Heaviside (1 = grounded, 0 = floating).
    rho_ratio : float
        rho_ice / rho_water.
    """
    rr = Constant(rho_ratio)
    return He * (b + H) + (Constant(1.0) - He) * (Constant(1.0) - rr) * H


def effective_pressure_ratio(H, b, rho_I=917.0, rho_W=1024.0, g=9.81, floor=0.01):
    """Effective pressure ratio phi_eff for sliding coefficient.

    phi_eff = max(1 - rho_W * g * max(0, -b) / (rho_I * g * H), floor)

    On floating ice phi_eff → floor (small), making K large → free slip.
    On grounded ice phi_eff → 1 (normal friction).

    Parameters
    ----------
    H : Function
        Ice thickness (should be > 0 everywhere).
    b : Function
        Bed elevation.
    floor : float
        Minimum phi_eff (prevents division by zero).
    """
    Q = H.function_space()
    rI = Constant(rho_I)
    rW = Constant(rho_W)
    grav = Constant(g)
    H_safe = max_value(H, Constant(1.0))  # protect against H=0
    return Function(Q).interpolate(
        max_value(
            Constant(1.0)
            - rW * grav * max_value(Constant(0.0), -b) / (rI * grav * H_safe),
            Constant(floor),
        )
    )


def sliding_coefficient(u_c, tau_c, phi_eff, n):
    """Weertman sliding coefficient with effective pressure.

    K = u_c / (phi_eff * tau_c)^n

    In icepack2's dual form: large K → free slip (floating), small K → friction.

    Parameters
    ----------
    u_c : Constant
        Reference velocity scale.
    tau_c : Constant
        Reference stress scale.
    phi_eff : Function
        Effective pressure ratio (0.01 on floating, 1 on grounded).
    n : Constant
        Glen/sliding exponent.
    """
    return u_c / (phi_eff * tau_c) ** n


def grounding_indicator(
    H, b, Q, method="smoothed_haf", kH=1.0, rho_I=917.0, rho_W=1024.0
):
    """Compute a grounding indicator for GL flux integration.

    Parameters
    ----------
    method : str
        'smoothed_haf' - smoothed height-above-flotation with variational
                         smoothing then sharp threshold (original approach)
        'heaviside' - smooth Heaviside He (differentiable, better for adjoint)
    """
    from firedrake import inner, grad, derivative, dx, solve

    rr = Constant(rho_I / rho_W)
    s = max_value(b + H, (Constant(1.0) - rr) * H)
    s_float = b + Constant(rho_W / rho_I) * max_value(-b, Constant(0.0))

    if method == "heaviside":
        haf = height_above_flotation(H, b, rho_I, rho_W)
        He = smooth_heaviside(haf, kH)
        return Function(Q).interpolate(He)

    elif method == "smoothed_haf":
        haf = Function(Q).interpolate(s - s_float)
        haf0 = haf.copy(deepcopy=True)
        solve(
            derivative(
                0.5
                * (
                    (haf - haf0) ** 2
                    + Constant(100.0) ** 2 * inner(grad(haf), grad(haf))
                )
                * dx,
                haf,
            )
            == 0,
            haf,
        )
        return Function(Q).interpolate(
            conditional(haf > 0, Constant(1.0), Constant(0.0))
        )
