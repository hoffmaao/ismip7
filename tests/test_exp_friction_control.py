r"""The exp friction control (ISMIP7_FRICTION_CONTROL=exp): C = C_ref exp(alpha)
with one scalar reference and a zero-mean prior on alpha. The residual the
inversion assembles (theta = alpha on a constant C_w0 = C_ref) is the residual
the log control assembles with theta = ln(C_ref exp(alpha) / C_w0) on the
anchor, so a MAP on either describes one friction; the start reproduces the
anchor exactly; and a step in alpha is a multiplicative change of C, the
same factor on a weak bed as on a stiff one (the reason for the control).

Serial, the slab of test_sqrt_friction_control, no data files.
"""
# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
from icepack2_tools.dual_friction import build_rc_residual     # noqa: E402

import numpy as np                                              # noqa: E402

from firedrake import (                                         # noqa: E402
    Constant, FiniteElement, Function, FunctionSpace, RectangleMesh,
    SpatialCoordinate, TensorFunctionSpace, VectorFunctionSpace, as_vector,
    assemble, conditional, exp, ln, max_value,
)

RHO_I, RHO_W = 917.0, 1024.0


def _slab():
    mesh = RectangleMesh(30, 10, 30e3, 10e3)
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Z = V * TensorFunctionSpace(mesh, dg0, symmetry=True) * VectorFunctionSpace(mesh, dg0)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, y = SpatialCoordinate(mesh)
    H = Function(Q0).interpolate(conditional(x < 20e3, 1500.0 - 0.03 * x, 400.0))
    b = Function(Q0).interpolate(conditional(x < 20e3, -200.0 - 0.01 * x, -1500.0))
    s = Function(Q0).interpolate(max_value(b + H, (1.0 - RHO_I / RHO_W) * H))
    # an anchor spanning two orders of magnitude: a weak bed along y = 0
    C_w0 = Function(Q, name="C_w0").interpolate(
        0.001 + 0.08 * (y / 10e3) ** 2 * (1.0 - x / 40e3))
    A_prior = Function(Q).interpolate(Constant(20.0) + 5.0 * x / 30e3)
    z = Function(Z)
    z.subfunctions[0].interpolate(as_vector([200.0 + 0.02 * x, 5.0 * y / 10e3]))
    z.subfunctions[1].interpolate(Constant(((0.05, 0.01), (0.01, 0.02))))
    z.subfunctions[2].interpolate(as_vector([0.03, 0.005]))
    return dict(mesh=mesh, Q=Q, Z=Z, H=H, b=b, s=s, C_w0=C_w0, A_prior=A_prior, z=z,
                n_flow=Constant(3.0), m_slide=Constant(3.0))


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


def test_the_start_reproduces_the_anchor_and_the_residual_is_the_log_controls():
    f = _slab()
    c_ref = 0.02
    # the inversion's start: alpha = ln(C_w0 / C_ref)
    alpha = Function(f["Q"]).interpolate(ln(max_value(f["C_w0"], Constant(1e-12)) / Constant(c_ref)))
    C_back = Function(f["Q"]).interpolate(Constant(c_ref) * exp(alpha)).dat.data_ro
    assert np.allclose(C_back, f["C_w0"].dat.data_ro, rtol=1e-12, atol=0.0)
    # the exp residual (C = C_ref exp(alpha) outright, theta = 0, as the
    # inversion assembles it) against the log residual on an anchor with
    # the deviation that describes the same friction. Uniform fields, as in
    # the sqrt test: the log path gates its deviation by the grounded
    # indicator, so on a varying anchor the two differ across the grounding
    # band by construction (the exp control has no anchor to fall back to).
    C_uniform = Function(f["Q"]).interpolate(Constant(0.05))
    alpha_u = Function(f["Q"]).interpolate(Constant(np.log(0.05 * 1.3 / c_ref)))
    r_exp = _residual(f, Constant(0.0), Constant(c_ref) * exp(alpha_u))
    dev = Function(f["Q"]).interpolate(ln(Constant(c_ref) * exp(alpha_u) / C_uniform))
    r_log = _residual(f, dev, C_uniform)
    scale = np.max(np.abs(r_log))
    assert scale > 0.0
    assert np.allclose(r_exp, r_log, rtol=0.0, atol=1e-9 * scale)


def test_a_step_in_alpha_is_the_same_factor_on_weak_and_stiff_beds():
    f = _slab()
    c_ref = 0.02
    alpha = Function(f["Q"]).interpolate(ln(max_value(f["C_w0"], Constant(1e-12)) / Constant(c_ref)))
    C0 = Function(f["Q"]).interpolate(Constant(c_ref) * exp(alpha)).dat.data_ro.copy()
    alpha.dat.data[:] -= 1.0
    C1 = Function(f["Q"]).interpolate(Constant(c_ref) * exp(alpha)).dat.data_ro
    ratio = C1 / C0
    assert np.allclose(ratio, np.exp(-1.0), rtol=1e-12)
    # under the sqrt control the same unit step is additive in sqrt(C): it
    # zeroes a weak bed (sqrt(0.001) = 0.03 < 1) and barely moves a stiff one
    a_sq = np.sqrt(C0)
    weak, stiff = a_sq < 0.05, a_sq > 0.25
    assert weak.any() and stiff.any()
    assert np.all(np.maximum(a_sq[weak] - 0.03, 0.0) ** 2 / C0[weak] < 0.5)
    assert np.all((a_sq[stiff] - 0.03) ** 2 / C0[stiff] > 0.75)
