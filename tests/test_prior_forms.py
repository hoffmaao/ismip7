r"""The two prior precisions in icepack2_tools.prior.

`laplacian` uses A itself; `bilaplacian` uses A M^-1 A, the
squared-operator ("LM^-1L") prior of Villa et al. (2021). These
check the second against dense linear algebra, because its energy is reached
through a mass solve rather than a single form and a sign or a transpose there
is invisible in the output of a run.
"""
import math

import numpy as np
import pytest

firedrake = pytest.importorskip("firedrake")
# Before any mesh: tlm_adjoint patches the function spaces built after it
# loads, and BilaplacianAuxSolver is a tlm_adjoint equation.
tlm = pytest.importorskip("tlm_adjoint.firedrake")
from firedrake import (  # noqa: E402
    Function, FunctionSpace, TestFunction, TrialFunction, UnitSquareMesh,
    assemble, dx, inner, solve,
)

from icepack2_tools.prior import (  # noqa: E402
    BilaplacianAuxSolver,
    bilaplacian_aux_residual, bilaplacian_coeffs, bilaplacian_energy_form,
    prior_bilinear_form, prior_operator_coeffs, prior_operator_form,
    regularization_form, regularization_gradient_form,
)

DELTA, GAMMA = 3.0, 0.7


@pytest.fixture
def space():
    return FunctionSpace(UnitSquareMesh(10, 10), "CG", 1)


@pytest.fixture
def theta(space):
    f = Function(space)
    f.dat.data[:] = np.random.default_rng(0).standard_normal(
        f.dat.data.shape[0]
    )
    return f


def _dense(space, delta, gamma):
    u, v = TrialFunction(space), TestFunction(space)
    A = np.array(assemble(prior_operator_form(u, v, delta, gamma)).M.handle[:, :])
    M = np.array(assemble(inner(u, v) * dx).M.handle[:, :])
    return A, M


def test_the_gradient_form_is_the_bilinear_form_at_the_control(theta, space):
    """regularization_gradient_form must stay A theta and not drift into its
    own copy of the operator: the metric and the UQ are built from the other
    one."""
    v = TestFunction(space)
    a = assemble(regularization_gradient_form(theta, v, 1e4, 1.0, L_reg=0.1))
    b = assemble(prior_bilinear_form(theta, v, 1e4, 1.0, L_reg=0.1))
    assert np.allclose(a.dat.data_ro, b.dat.data_ro, rtol=0, atol=0)


def test_the_laplacian_energy_is_half_theta_A_theta(theta, space):
    delta, gamma = prior_operator_coeffs(1e4, 2.0, 0.1)
    A, _ = _dense(space, delta, gamma)
    x = theta.dat.data_ro
    got = float(assemble(regularization_form(theta, 1e4, 2.0, L_reg=0.1)))
    assert got == pytest.approx(0.5 * x @ A @ x, rel=1e-12)


def test_the_bilaplacian_energy_is_half_theta_A_Minv_A_theta(theta, space):
    """The squared precision, against dense algebra."""
    A, M = _dense(space, DELTA, GAMMA)
    B = A @ np.linalg.solve(M, A)
    aux = Function(space)
    solve(
        bilaplacian_aux_residual(theta, aux, TestFunction(space), DELTA, GAMMA)
        == 0,
        aux,
    )
    got = float(assemble(bilaplacian_energy_form(aux)))
    x = theta.dat.data_ro
    assert got == pytest.approx(0.5 * x @ B @ x, rel=1e-10)


def test_the_bilaplacian_gradient_is_A_times_the_solved_field(theta, space):
    """dR/dtheta = A M^-1 A theta = A f, which is why the inversion reuses the
    energy's f rather than solving twice."""
    A, M = _dense(space, DELTA, GAMMA)
    v = TestFunction(space)
    aux = Function(space)
    solve(bilaplacian_aux_residual(theta, aux, v, DELTA, GAMMA) == 0, aux)
    got = assemble(prior_operator_form(aux, v, DELTA, GAMMA)).dat.data_ro
    want = A @ np.linalg.solve(M, A @ theta.dat.data_ro)
    assert np.allclose(got, want, rtol=1e-9, atol=1e-12)


def test_the_factored_aux_solve_is_the_residual_solve(theta, space):
    """BilaplacianAuxSolver, which the inversion uses off the tape, against
    the residual solve it replaces and dense algebra: two controls and two
    (delta, gamma) pairs through one factor, so a stale right-hand side or a
    factor of the wrong operator cannot pass on the first call alone."""
    v = TestFunction(space)
    solver = BilaplacianAuxSolver(
        space, form_compiler_parameters={"quadrature_degree": 4})
    factor = None
    for scale, (delta, gamma) in ((1.0, (DELTA, GAMMA)), (-2.5, (1.3, 0.2))):
        A, M = _dense(space, delta, gamma)
        ctrl = Function(space).assign(scale * theta)
        ref = Function(space)
        solve(bilaplacian_aux_residual(ctrl, ref, v, delta, gamma) == 0, ref)
        got = solver(ctrl, Function(space), delta, gamma)
        assert np.allclose(got.dat.data_ro, ref.dat.data_ro, rtol=1e-10, atol=1e-12)
        want = np.linalg.solve(M, A @ ctrl.dat.data_ro)
        assert np.allclose(got.dat.data_ro, want, rtol=1e-9, atol=1e-12)
        factor = factor or solver._mats
        assert solver._mats is factor


def test_the_taped_factored_solve_has_the_bilaplacian_gradient(theta, space):
    """On the tape (the TAO path's objective) the factored solve must give
    dR/dtheta = A M^-1 A theta, as the residual-form EquationSolver it
    replaces does, and pass tlm_adjoint's Taylor tests."""
    A, M = _dense(space, DELTA, GAMMA)
    solver = BilaplacianAuxSolver(space)

    def energy_factored(ctrl):
        aux = Function(space)
        solver(ctrl, aux, DELTA, GAMMA)
        J = tlm.Functional(name="R")
        J.assign(bilaplacian_energy_form(aux))
        return J

    def energy_residual(ctrl):
        aux = Function(space)
        tlm.EquationSolver(
            bilaplacian_aux_residual(ctrl, aux, TestFunction(space), DELTA, GAMMA)
            == 0, aux).solve()
        J = tlm.Functional(name="R")
        J.assign(bilaplacian_energy_form(aux))
        return J

    grads = {}
    for name, fwd in (("factored", energy_factored), ("residual", energy_residual)):
        tlm.reset_manager()
        tlm.start_manager()
        J = fwd(theta)
        tlm.stop_manager()
        grads[name] = (float(J), tlm.compute_gradient(J, theta))
    J_val, dJ = grads["factored"]
    x = theta.dat.data_ro
    want = A @ np.linalg.solve(M, A @ x)
    assert J_val == pytest.approx(0.5 * x @ want, rel=1e-10)
    assert np.allclose(dJ.dat.data_ro, want, rtol=1e-9, atol=1e-12)
    assert np.allclose(dJ.dat.data_ro, grads["residual"][1].dat.data_ro,
                       rtol=1e-9, atol=1e-12)
    # A wrong adjoint gives order 1; the energy is quadratic, so order 2.
    assert tlm.taylor_test(energy_factored, theta, J_val=J_val, dJ=dJ) > 1.99
    assert tlm.taylor_test_tlm(energy_factored, theta, tlm_order=1) > 1.99
    tlm.reset_manager()


def test_the_covariance_action_inverts_the_precision(space):
    """A^-1 M A^-1 is the inverse of A M^-1 A: the metric and the prior are
    the same object, not two operators that merely look alike."""
    A, M = _dense(space, DELTA, GAMMA)
    B = A @ np.linalg.solve(M, A)
    g = np.random.default_rng(1).standard_normal(A.shape[0])
    cov = np.linalg.solve(A, M @ np.linalg.solve(A, g))
    assert np.allclose(B @ cov, g, rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize("sigma,rho", [(0.3, 7500.0), (1.0, 2500.0), (0.05, 500.0)])
def test_sigma_and_rho_round_trip_through_the_villa_2021_relations(sigma, rho):
    """sigma^2 = 1/(4 pi gamma delta) and rho = sqrt(8 gamma/delta) at nu=1."""
    delta, gamma = bilaplacian_coeffs(sigma, rho)
    assert 1.0 / math.sqrt(4 * math.pi * gamma * delta) == pytest.approx(sigma, rel=1e-12)
    assert math.sqrt(8 * gamma / delta) == pytest.approx(rho, rel=1e-12)


@pytest.mark.parametrize("sigma,rho", [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0)])
def test_a_degenerate_prior_is_refused(sigma, rho):
    with pytest.raises(ValueError):
        bilaplacian_coeffs(sigma, rho)
