"""prior.py - Whittle-Matern prior for the ISMIP7 inversion.

Ported from mismip_time-dependent-da/prior.py (Recinos et al. 2023
convention). For a control field ``theta`` (= log(param / param_prior), a
log-deviation from a PHYSICAL prior mean) with strength ``gamma`` and
correlation length ``L_reg``, the prior energy added to the inversion cost is

    R(theta; gamma) = 0.5/area * gamma * \\int [ theta^2 + L_reg^2 |grad theta|^2 ] dx

Its Hessian is the prior precision A = delta_eff*M + gamma_eff*K with
delta_eff = gamma/area (mass) and gamma_eff = gamma*L_reg^2/area (stiffness),
so the correlation length sqrt(gamma_eff/delta_eff) = L_reg is fixed and gamma
is the single per-control strength knob.

WHY THIS FIXES THE ISMIP7 n=3 BLOW-UP: the old ISMIP7 regularization was
PURE SMOOTHNESS (gamma*L^2*|grad theta|^2 only, no mass term), so a large
smooth control field was unpenalized -- and with the fluidity control on a
CONSTANT baseline A0, phi had to carry all the spatial fluidity structure and
blew up to +-36 at n=3. The fix of Recinos et al. (2023) is two-fold and does NOT penalize
parameter amplitude: (1) put the control on a PHYSICAL prior mean (thermo
fluidity A_prior, balance friction C_w0), so the amplitude lives in the mean
and theta is a small deviation; (2) use this proper prior whose mass term just
removes the null space, at a physical variance (gamma ~ 1e4 -> prior std
~0.2-0.3, the physical log-deviation scale). The regularization then constrains
the DEVIATION, not the amplitude.
"""

L_REG = 7500.0  # correlation length (m); the one definition used everywhere


def regularization_form(theta, gamma, area, L_reg=L_REG):
    """Whittle-Matern prior energy 0.5/area * gamma * (theta^2 + L_reg^2
    |grad theta|^2) dx, as a UFL form. Add to the inversion cost per control."""
    from firedrake import dx, inner, grad
    coef = 0.5 * float(gamma) / float(area)
    return coef * (theta ** 2
                   + float(L_reg) ** 2 * inner(grad(theta), grad(theta))) * dx


def prior_operator_coeffs(gamma, area, L_reg=L_REG):
    """``(delta_eff, gamma_eff)`` of the prior precision ``A = delta_eff*M +
    gamma_eff*K``: exactly the coefficients of the Hessian of
    :func:`regularization_form`."""
    gamma = float(gamma)
    area = float(area)
    return gamma / area, gamma * float(L_reg) ** 2 / area


def prior_bilinear_form(trial, test, gamma, area, L_reg=L_REG):
    """Prior precision ``A`` as a bilinear form.

    The Hessian of :func:`regularization_form`, so the regularization the
    inversion pays, the metric it descends in, and the operator the UQ
    eigendecomposes against are one definition rather than three copies that
    can drift (the reason the MISMIP TDDA pipeline keeps them in one module).
    """
    from firedrake import dx, inner, grad
    delta_eff, gamma_eff = prior_operator_coeffs(gamma, area, L_reg)
    return (delta_eff * inner(trial, test)
            + gamma_eff * inner(grad(trial), grad(test))) * dx


def regularization_gradient_form(theta, test, gamma, area, L_reg=L_REG):
    """Gateaux derivative of regularization_form wrt theta (the assembled
    cost-gradient contribution): ``A theta`` tested against ``test``, which is
    :func:`prior_bilinear_form` evaluated at the current control."""
    return prior_bilinear_form(theta, test, gamma, area, L_reg=L_reg)


# ---------------------------------------------------------------------------
# Bi-Laplacian (squared) prior -- the precision of Villa et al. (2021)
# ---------------------------------------------------------------------------
# Everything above uses A ITSELF as the prior precision. Villa et al. (2021)
# use the SQUARE ("LM^-1L") instead,
#
#     B = A M^-1 A,        A = delta*M + gamma*K
#
# and that is not a cosmetic difference. The Whittle-Matern SPDE
# (delta - gamma*Lap)^(alpha/2) u = W gives a field in L^2 only for
# alpha > d/2; in d = 2, A alone is alpha = 1 and the "field" is a
# distribution, not a function -- its pointwise variance does not exist and
# refining the mesh does not converge to anything. Squaring gives alpha = 2,
# nu = alpha - d/2 = 1, and a genuine function-valued Matern field. This is
# why both codes square, and why both quote closed-form marginal variance and
# correlation length, which the un-squared operator has none of.
#
# The energy needs a mass solve, so unlike regularization_form it is not one
# UFL form: R(theta) = 0.5 * theta' A M^-1 A theta is evaluated by solving
# M f = A theta and then integrating 0.5*f^2 (exactly the
# norm_sq applied to its solved field). The caller owns the solve so that the
# inversion can put it on the adjoint tape and the UQ can reuse the operator.

def bilaplacian_coeffs(sigma, rho):
    r"""``(delta, gamma)`` of ``A = delta*M + gamma*K`` for a 2-D Matern field
    of marginal standard deviation ``sigma`` and correlation length ``rho`` (m)
    under the SQUARED precision ``A M^-1 A``.

    the relations of Villa et al. (2021) at ``nu = alpha - d/2 = 1``:

        sigma^2 = 1 / (4*pi*gamma*delta),      rho = sqrt(8*gamma/delta)

    inverted here. Unlike the single-knob ``gamma``/``L_reg`` pair of the
    un-squared form, both numbers are physical: ``sigma`` is the log-deviation
    scale the control is expected to carry and ``rho`` the distance over which
    it decorrelates, so a prior can be SET rather than tuned.
    """
    import math
    sigma = float(sigma)
    rho = float(rho)
    if sigma <= 0.0 or rho <= 0.0:
        raise ValueError(f"sigma and rho must be positive, got {sigma}, {rho}")
    delta = math.sqrt(2.0 / math.pi) / (sigma * rho)
    gamma = delta * rho ** 2 / 8.0
    return delta, gamma


def prior_operator_form(trial, test, delta, gamma):
    r"""``A = delta*M + gamma*K`` as a bilinear form, in the (delta, gamma)
    parameterisation Villa et al. (2021) use directly.

    :func:`prior_bilinear_form` is the same operator reached through this
    repository's (gamma, area, L_reg) knobs; both exist so the two conventions
    never drift into two different operators.
    """
    from firedrake import dx, inner, grad
    return (float(delta) * inner(trial, test)
            + float(gamma) * inner(grad(trial), grad(test))) * dx


def bilaplacian_aux_residual(theta, aux, test, delta, gamma):
    r"""Residual of the mass solve ``M f = A theta`` defining ``f = M^-1 A theta``.

    Solve this for ``aux`` (inside the tape, so the adjoint carries it), then
    the prior energy is :func:`bilaplacian_energy_form` of the result.
    """
    from firedrake import dx, inner
    return (inner(aux, test) * dx
            - prior_operator_form(theta, test, delta, gamma))


def bilaplacian_energy_form(aux):
    r"""``0.5 * \int f^2 dx`` with ``f = M^-1 A theta``, i.e.
    ``0.5 * theta' A M^-1 A theta`` -- the squared norm of the solved
    field."""
    from firedrake import dx, inner
    return 0.5 * inner(aux, aux) * dx
