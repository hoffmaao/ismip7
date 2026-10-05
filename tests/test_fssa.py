"""Free-surface stabilization term (icepack2_tools/fssa.py)."""
import numpy as np
import firedrake as fd
from firedrake import (Constant, Function, FunctionSpace, VectorFunctionSpace,
                       TensorFunctionSpace, SpatialCoordinate, assemble, derivative,
                       as_vector)

from icepack2_tools.fssa import fssa_term
from icepack_tools.constants import gravity, ice_density, water_density


def _setup(bed):
    mesh = fd.UnitSquareMesh(6, 6)
    mesh.coordinates.dat.data[:] *= 1e4
    V = VectorFunctionSpace(mesh, "CG", 1)
    Z = V * TensorFunctionSpace(mesh, "DG", 0, symmetry=True) * V
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, y = SpatialCoordinate(mesh)
    z = Function(Z)
    z.subfunctions[0].interpolate(as_vector([100.0 + 0.02 * x, 0.01 * y]))
    H = Function(Q0).interpolate(Constant(800.0) + 0.001 * x)
    b = Function(Q0).interpolate(Constant(bed))
    return mesh, z, H, b


def _u_block(F, z):
    """The velocity block of the Jacobian, dense (one rank: the first V.dim() dofs)."""
    n = z.subfunctions[0].function_space().dim()
    J = assemble(derivative(F, z), mat_type="aij").petscmat
    return J[:n, :n]


def test_vanishes_at_the_reference_and_without_a_step():
    mesh, z, H, b = _setup(100.0)
    u_ref = Function(z.subfunctions[0].function_space()).assign(z.subfunctions[0])
    assert fssa_term(z, u_ref, None, H, b) == 0 and fssa_term(z, None, Constant(0.1), H, b) == 0
    r = assemble(fssa_term(z, u_ref, Constant(0.1), H, b)).dat.data_ro
    assert np.abs(np.concatenate([a.ravel() for a in r])).max() == 0.0


def test_symmetric_semidefinite_and_scaled_by_the_step():
    mesh, z, H, b = _setup(100.0)
    u_ref = Function(z.subfunctions[0].function_space())
    J1 = _u_block(fssa_term(z, u_ref, Constant(0.05), H, b), z)
    J2 = _u_block(fssa_term(z, u_ref, Constant(0.1), H, b), z)
    assert np.allclose(J1, J1.T) and np.allclose(2 * J1, J2)
    # the residual's convention is the negative of the energy form, so the
    # block is negative semi-definite, like the membrane term it joins
    assert np.linalg.eigvalsh(J1).max() < 1e-12 * np.abs(J1).max()


def test_floating_ice_carries_the_surface_ratio():
    """On grounded ice the surface moves with the thickness; afloat only
    1 - rho_i/rho_w of it does, and the stabilization follows."""
    _, zg, Hg, bg = _setup(100.0)
    _, zf, Hf, bf = _setup(-5000.0)
    Jg = _u_block(fssa_term(zg, Function(zg.subfunctions[0].function_space()), Constant(0.1), Hg, bg), zg)
    Jf = _u_block(fssa_term(zf, Function(zf.subfunctions[0].function_space()), Constant(0.1), Hf, bf), zf)
    assert np.allclose(Jf, (1.0 - ice_density / water_density) * Jg, rtol=1e-10)


def test_magnitude_is_rho_g_h2_grad_div():
    mesh, z, H, b = _setup(100.0)
    V = z.subfunctions[0].function_space()
    u_ref = Function(V)
    u, v = fd.TrialFunction(V), fd.TestFunction(V)
    expected = -assemble(Constant(0.1 * ice_density * gravity) * H ** 2 * fd.div(u) * fd.div(v) * fd.dx, mat_type="aij").petscmat[:, :]
    J = _u_block(fssa_term(z, u_ref, Constant(0.1), H, b), z)
    assert np.allclose(J, expected, rtol=1e-10)
