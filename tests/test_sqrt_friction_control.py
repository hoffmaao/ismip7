r"""The sqrt(C) friction control (ISMIP7_FRICTION_CONTROL=sqrt): the residual
the inversion assembles with ``C_w0 = alpha^2`` and a zero log deviation is
the residual the log control assembles with ``theta = ln(alpha^2 / C_w0)`` on
the anchor, so the two parameterisations describe one friction and a MAP on
either can be read by a consumer of the other. Also that the zero-mean prior
on alpha is the same bi-Laplacian operator, at the friction's own scale.

Serial, a 30 km x 10 km grounded slab with a floating tongue, no data files.
"""
# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
from icepack2_tools.dual_friction import build_rc_residual     # noqa: E402

import numpy as np                                              # noqa: E402
import pytest                                                   # noqa: E402

from firedrake import (                                         # noqa: E402
    Constant,
    FiniteElement,
    Function,
    FunctionSpace,
    RectangleMesh,
    SpatialCoordinate,
    TensorFunctionSpace,
    VectorFunctionSpace,
    as_vector,
    assemble,
    conditional,
    exp,
    ln,
    max_value,
    sqrt,
)

from icepack2_tools.prior import bilaplacian_coeffs             # noqa: E402

RHO_I, RHO_W = 917.0, 1024.0


@pytest.fixture(scope="module")
def slab():
    mesh = RectangleMesh(30, 10, 30e3, 10e3)
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Z = V * TensorFunctionSpace(mesh, dg0, symmetry=True) * VectorFunctionSpace(mesh, dg0)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, y = SpatialCoordinate(mesh)
    # grounded for x < 20 km, floating beyond
    H = Function(Q0).interpolate(conditional(x < 20e3, 1500.0 - 0.03 * x, 400.0))
    b = Function(Q0).interpolate(conditional(x < 20e3, -200.0 - 0.01 * x, -1500.0))
    s = Function(Q0).interpolate(
        max_value(b + H, (1.0 - RHO_I / RHO_W) * H))
    # a smooth anchor that varies along and across flow, on the control's
    # own space so ln(alpha^2 / C_w0) is well defined at every node
    C_w0 = Function(Q, name="C_w0").interpolate(
        0.05 + 0.03 * (1.0 + y / 10e3) * (1.0 - x / 40e3))
    A_prior = Function(Q).interpolate(Constant(20.0) + 5.0 * x / 30e3)
    z = Function(Z)
    z.subfunctions[0].interpolate(as_vector([200.0 + 0.02 * x, 5.0 * y / 10e3]))
    z.subfunctions[1].interpolate(Constant(((0.05, 0.01), (0.01, 0.02))))
    z.subfunctions[2].interpolate(as_vector([0.03, 0.005]))
    n_flow, m_slide = Constant(3.0), Constant(3.0)
    return dict(mesh=mesh, Q=Q, Q0=Q0, Z=Z, H=H, b=b, s=s, C_w0=C_w0,
                A_prior=A_prior, z=z, n_flow=n_flow, m_slide=m_slide)


def _residual(f, theta, C):
    F = build_rc_residual(
        f["z"], theta, Function(f["Q"]), H=f["H"], s=f["s"], b=f["b"], C_w0=C,
        A4_base=f["A_prior"], n_flow=f["n_flow"], n_flow_val=3.0,
        m_slide=f["m_slide"], tau_c=Constant(0.1), alpha=Constant(1e-2),
        H_ref=Constant(100.0), fric_law="budd", N_ref=None,
        nhat_floor=0.02, nhat_cap=3.0, alpha_gl=0.5, c_w0_floor=0.0,
        h_visc_floor=10.0, ocean_drag=0.0, k_lim=0.0,
    )
    with assemble(F).dat.vec_ro as v:
        return v.array.copy()


def test_alpha_squared_is_the_log_control_on_the_anchor(slab):
    f = slab
    # Uniform fields: the two parameterisations coincide pointwise, so the
    # residuals agree to roundoff. (With a varying alpha the nodal
    # interpolant of ln(alpha^2 / C_w0) no longer squares to alpha^2 inside
    # a cell, which is a representation difference, not a physics one.)
    C_uniform = Function(f["Q0"]).interpolate(Constant(0.05))
    alpha = Function(f["Q"], name="alpha").interpolate(Constant(np.sqrt(0.05) * 1.3))
    r_sqrt = _residual(f, Constant(0.0), alpha ** 2)
    theta = Function(f["Q"]).interpolate(ln(alpha ** 2 / C_uniform))
    r_log = _residual(f, theta, C_uniform)
    scale = np.max(np.abs(r_log))
    assert scale > 0.0
    assert np.allclose(r_sqrt, r_log, rtol=0.0, atol=1e-9 * scale)


def test_a_varying_alpha_agrees_with_the_log_control_to_second_order(slab):
    f = slab
    x, y = SpatialCoordinate(f["mesh"])
    alpha = Function(f["Q"], name="alpha").interpolate(
        sqrt(f["C_w0"]) * (1.2 + 0.3 * y / 10e3))
    r_sqrt = _residual(f, Constant(0.0), alpha ** 2)
    theta = Function(f["Q"]).interpolate(ln(alpha ** 2 / f["C_w0"]))
    r_log = _residual(f, theta, f["C_w0"])
    scale = np.max(np.abs(r_log))
    # 30% variation over 10 km on 1 km cells: the in-cell mismatch of
    # exp(interpolated log) against the squared interpolant is second order
    diff = np.max(np.abs(r_sqrt - r_log)) / scale
    assert diff < 1e-3, diff
    assert not np.allclose(r_sqrt, r_log, rtol=0.0, atol=1e-12 * scale)


def test_the_friction_it_describes_is_alpha_squared(slab):
    f = slab
    alpha = Function(f["Q"]).interpolate(sqrt(f["C_w0"]) * 1.5)
    r_a = _residual(f, Constant(0.0), alpha ** 2)
    # a different alpha gives a different residual: the control is live
    r_b = _residual(f, Constant(0.0), (alpha * 0.5) ** 2)
    assert not np.allclose(r_a, r_b, rtol=0.0, atol=1e-9 * np.max(np.abs(r_a)))
    # alpha^2 is non-negative whatever sign alpha takes, so a sign flip is
    # the same friction (Recinos's parameterisation needs no bound on alpha)
    r_c = _residual(f, Constant(0.0), (-alpha) ** 2)
    assert np.allclose(r_a, r_c, rtol=0.0, atol=1e-12 * np.max(np.abs(r_a)))


def test_the_prior_scale_is_the_friction_scale():
    # sigma in alpha units: the operator's coefficients follow the same
    # closed forms as the log control's, at the friction's own scale
    d1, g1 = bilaplacian_coeffs(0.127, 7500.0)
    d2, g2 = bilaplacian_coeffs(0.3, 7500.0)
    assert d1 / d2 == pytest.approx(0.3 / 0.127)
    assert g1 / g2 == pytest.approx(0.3 / 0.127)
    assert np.sqrt(8.0 * g1 / d1) == pytest.approx(7500.0)
