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

import functools

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
# why the squared form is used, and why it has closed-form marginal variance
# and correlation length, which the un-squared operator has none of.
#
# The energy needs a mass solve, so unlike regularization_form it is not one
# UFL form: R(theta) = 0.5 * theta' A M^-1 A theta is evaluated by solving
# M f = A theta and then integrating 0.5*f^2 (the squared L2
# norm of the solved field). The caller owns the solve so that the
# inversion can put it on the adjoint tape and the UQ can reuse the operator.

def bilaplacian_coeffs(sigma, rho):
    r"""``(delta, gamma)`` of ``A = delta*M + gamma*K`` for a 2-D Matern field
    of marginal standard deviation ``sigma`` and correlation length ``rho`` (m)
    under the SQUARED precision ``A M^-1 A``.

    The relations of Villa et al. (2021) at ``nu = alpha - d/2 = 1``:

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


# MUMPS Cholesky, the options the inversion's prior-metric solvers use. M is
# symmetric positive definite on any mesh.
MASS_CHOLESKY = {"ksp_type": "preonly", "pc_type": "cholesky",
                 "pc_factor_mat_solver_type": "mumps"}


def _update(vec, other, method):
    if method == "assign":
        other.copy(vec)
    elif method == "add":
        vec.axpy(1.0, other)
    elif method == "sub":
        vec.axpy(-1.0, other)
    else:
        raise ValueError(f"unknown method {method!r}")


class BilaplacianAuxSolver:
    r"""``f = M^-1 A theta`` for the bi-Laplacian energy, with ``M`` factored
    on the first call and every later call a back-substitution.

    ``M`` and the stiffness ``K`` depend on the mesh alone, so one factor and
    one pair of matrices serve both controls (``A = delta*M + gamma*K``) and
    every evaluation. Each call is a tlm_adjoint ``LinearEquation``: recorded
    while a manager annotates (the TAO path's objective, whose adjoint
    back-substitutes with the same factor, ``M`` being symmetric), an ordinary
    solve when none does.

    Solving :func:`bilaplacian_aux_residual` instead takes Firedrake's default
    Newton solve, a fresh MUMPS LU of ``M`` at every call (and a second one in
    the adjoint), whose analysis MUMPS runs on one rank by default
    (ICNTL(28)=1). Its cost grows with the vertex count and not with the rank
    count. For both controls at 1.2 million vertices that was 8.5 to 9.4 s an
    evaluation on 1 to 8 ranks of a workstation, against 0.12 to 0.23 s for
    the back-substitution after a 4.2 s factorisation
    (antarctica/scripts/probe_eval_overhead.py; antarctica/README.md,
    "Inversion solver"). tlm_adjoint's own caches cannot hold the factor:
    they are cleared before every evaluation, by the driver's forward and by
    TAOSolver's ReducedFunctional.

    The answer is the residual solve's to rounding: the same ``M`` and ``A``
    at the same quadrature, and an exact factorisation. tlm_adjoint must be
    imported before the mesh is built, as for any tlm_adjoint equation.
    """

    def __init__(self, space, *, form_compiler_parameters=None):
        self._space = space
        self._fcp = dict(form_compiler_parameters or {})
        self._mats = None

    def _setup(self):
        from firedrake import (LinearSolver, TestFunction, TrialFunction,
                               assemble, dx, grad, inner)
        if self._mats is None:
            u, v = TrialFunction(self._space), TestFunction(self._space)
            M = assemble(inner(u, v) * dx, form_compiler_parameters=self._fcp)
            K = assemble(inner(grad(u), grad(v)) * dx,
                         form_compiler_parameters=self._fcp)
            self._mats = (M.petscmat, K.petscmat,
                          LinearSolver(M, solver_parameters=MASS_CHOLESKY))
        return self._mats

    def __call__(self, theta, aux, delta, gamma):
        r"""Solve ``M aux = A theta`` into ``aux``; returns ``aux``."""
        from tlm_adjoint.firedrake import LinearEquation
        M, K, solver = self._setup()
        LinearEquation(
            aux, _PriorActionRHS(theta, M, K, float(delta), float(gamma)),
            A=_FactoredMass(self._space, M, solver),
        ).solve()
        return aux


@functools.lru_cache(maxsize=None)
def _matrix_classes():
    # Firedrake and tlm_adjoint are imported lazily everywhere in this
    # module, so a caller that needs only bilaplacian_coeffs runs without
    # them.
    from firedrake import Function
    from tlm_adjoint.firedrake import Matrix, RHS

    class FactoredMass(Matrix):
        r"""``M``, solved by back-substitution with a factor held outside
        tlm_adjoint's caches. Symmetric, so the adjoint solve is the forward
        one."""

        def __init__(self, space, M, solver):
            super().__init__(nl_deps=[], ic=False, adj_ic=False)
            self._space, self._M, self._solver = space, M, solver

        def forward_action(self, nl_deps, x, b, *, method="assign"):
            with x.dat.vec_ro as xv, b.dat.vec as bv:
                y = bv.duplicate()
                self._M.mult(xv, y)
                _update(bv, y, method)

        def adjoint_action(self, nl_deps, adj_x, b, b_index=0, *, method="assign"):
            if b_index != 0:
                raise ValueError("unexpected b_index")
            self.forward_action(nl_deps, adj_x, b, method=method)

        def forward_solve(self, x, nl_deps, b):
            self._solver.solve(x, b)

        def adjoint_solve(self, adj_x, nl_deps, b):
            if adj_x is None:
                adj_x = Function(self._space)
            self._solver.solve(adj_x, b)
            return adj_x

        def tangent_linear_rhs(self, tlm_map, x):
            return None

    class PriorActionRHS(RHS):
        r"""``A theta = delta*M theta + gamma*K theta``, linear in ``theta``."""

        def __init__(self, theta, M, K, delta, gamma):
            super().__init__([theta], nl_deps=[])
            self._M, self._K, self._dg = M, K, (delta, gamma)

        def _act(self, x, b, method):
            with x.dat.vec_ro as xv, b.dat.vec as bv:
                y = bv.duplicate()
                t = bv.duplicate()
                self._M.mult(xv, y)
                self._K.mult(xv, t)
                y.scale(self._dg[0])
                y.axpy(self._dg[1], t)
                _update(bv, y, method)

        def add_forward(self, b, deps):
            (theta,) = deps
            self._act(theta, b, "add")

        def subtract_adjoint_derivative_action(self, nl_deps, dep_index, adj_x, b):
            if dep_index != 0:
                raise ValueError("unexpected dep_index")
            # A is symmetric: A^* adj_x = A adj_x
            self._act(adj_x, b, "sub")

        def tangent_linear_rhs(self, tlm_map):
            tau = tlm_map[self.dependencies()[0]]
            if tau is None:
                return None
            return PriorActionRHS(tau, self._M, self._K, *self._dg)

    return FactoredMass, PriorActionRHS


def _FactoredMass(space, M, solver):
    return _matrix_classes()[0](space, M, solver)


def _PriorActionRHS(theta, M, K, delta, gamma):
    return _matrix_classes()[1](theta, M, K, delta, gamma)
