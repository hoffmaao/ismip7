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


def _velocity_residual(form):
    """The velocity part of an assembled residual, on a mixed space or on V."""
    data = assemble(form).dat.data_ro
    return np.asarray(data[0] if isinstance(data, tuple) else data).ravel()


def test_the_tendency_adds_a_load_and_leaves_the_jacobian():
    """``T(u_ref)`` only moves the right-hand side: the step reference is as
    stable as the start reference, and its load is tau rho g gamma h T div(v)
    (gamma = 1 on this grounded slab)."""
    mesh, z, H, b = _setup(100.0)
    V = z.subfunctions[0].function_space()
    u_ref = Function(V)
    x, y = SpatialCoordinate(mesh)
    T = Function(FunctionSpace(mesh, "DG", 0)).interpolate(3.0 - 2e-4 * x + 1e-4 * y)
    tau = Constant(0.05)
    with_t = fssa_term(z, u_ref, tau, H, b, tendency=T)
    without = fssa_term(z, u_ref, tau, H, b)
    assert np.allclose(_u_block(with_t, z), _u_block(without, z), rtol=0, atol=0)
    v = fd.TestFunction(V)
    load = _velocity_residual(Constant(0.05 * ice_density * gravity) * H * T * fd.div(v) * fd.dx)
    diff = _velocity_residual(with_t) - _velocity_residual(without)
    assert np.allclose(diff, load, rtol=1e-12, atol=1e-12 * np.abs(load).max())


def _upwind(mesh, u, h):
    n = fd.FacetNormal(mesh)
    un = fd.dot(u, n)
    unp = (un + abs(un)) / 2
    phi = fd.TestFunction(h.function_space())
    return ((unp("+") * h("+") - unp("-") * h("-")) * fd.jump(phi) * fd.dS
            + unp * h * phi * fd.ds)


def _advance(mesh, u, h_old, src, dt):
    """One implicit upwind DG0 advance, the forward's transport form."""
    Q0 = h_old.function_space()
    h_trial = fd.TrialFunction(Q0)
    phi = fd.TestFunction(Q0)
    n = fd.FacetNormal(mesh)
    un = fd.dot(u, n)
    unp = (un + abs(un)) / 2
    F = ((h_trial - h_old) / Constant(dt) * phi * fd.dx
         + (unp("+") * h_trial("+") - unp("-") * h_trial("-")) * fd.jump(phi) * fd.dS
         + unp * h_trial * phi * fd.ds - src * phi * fd.dx)
    h_new = Function(Q0)
    fd.solve(fd.lhs(F) == fd.rhs(F), h_new,
             solver_parameters={"ksp_type": "preonly", "pc_type": "lu"})
    return h_new


def _cell_area(Q0):
    return assemble(fd.TestFunction(Q0) * fd.dx).dat.data_ro


def test_the_realized_tendency_is_the_source_less_the_new_divergence():
    """The step reference reads T off the advance as (h_new - h_old) / dt; the
    implicit upwind transport makes that the source less the flux divergence
    of the new thickness under the velocity it used."""
    mesh, z, H, b = _setup(100.0)
    u = z.subfunctions[0]
    Q0 = H.function_space()
    x, y = SpatialCoordinate(mesh)
    src = Function(Q0).interpolate(0.3 - 5e-5 * x + 2e-5 * y)
    dt = 0.05
    h_new = _advance(mesh, u, H, src, dt)
    T = (h_new.dat.data_ro - H.dat.data_ro) / dt
    div_new = assemble(_upwind(mesh, u, h_new)).dat.data_ro / _cell_area(Q0)
    assert np.allclose(T, src.dat.data_ro - div_new, rtol=0, atol=1e-8 * np.abs(T).max())
    assert np.abs(T).max() > 1.0


def test_the_step_reference_vanishes_at_a_steady_state():
    """Fed the source that balances the flux divergence, the advance leaves
    the thickness where it is, T = 0, and with u = u_ref the step-referenced
    term is zero: no load is left at a steady state, whatever its forcing."""
    mesh, z, H, b = _setup(100.0)
    u = z.subfunctions[0]
    Q0 = H.function_space()
    src = Function(Q0)
    src.dat.data[:] = assemble(_upwind(mesh, u, H)).dat.data_ro / _cell_area(Q0)
    h_new = _advance(mesh, u, H, src, 0.05)
    T = Function(Q0)
    T.dat.data[:] = (h_new.dat.data_ro - H.dat.data_ro) / 0.05
    assert np.abs(T.dat.data_ro).max() < 1e-8 * np.abs(src.dat.data_ro).max()
    u_ref = Function(u.function_space()).assign(u)
    r = _velocity_residual(fssa_term(z, u_ref, Constant(0.05), H, b, tendency=T))
    # the start reference measured from a velocity of another divergence
    # leaves a load at the same state; the step reference removes it
    x, y = SpatialCoordinate(mesh)
    u_start = Function(u.function_space()).interpolate(
        as_vector([100.0 + 0.05 * x, 0.01 * y]))
    r_start = _velocity_residual(fssa_term(z, u_start, Constant(0.05), H, b))
    assert np.abs(r).max() < 1e-6 * np.abs(r_start).max()


def test_auto_follows_the_apparent_mass_balance():
    from icepack2_tools.fssa import resolve_reference
    assert resolve_reference("auto", True) == "start"
    assert resolve_reference("auto", False) == "step"
    for ref in ("start", "step"):
        assert resolve_reference(ref, True) == ref == resolve_reference(ref, False)
    import pytest
    with pytest.raises(ValueError, match="ISMIP7_FSSA_REFERENCE"):
        resolve_reference("previous", False)


def test_the_knob_defaults_to_auto(monkeypatch):
    from icepack2_tools import solverconfig
    monkeypatch.delenv("ISMIP7_FSSA_REFERENCE", raising=False)
    assert solverconfig.fssa_reference() == "auto"
    monkeypatch.setenv("ISMIP7_FSSA_REFERENCE", " Step ")
    assert solverconfig.fssa_reference() == "step"
    monkeypatch.setenv("ISMIP7_FSSA_REFERENCE", "")
    assert solverconfig.fssa_reference() == "auto"
    assert solverconfig.effective_solver_env()["ISMIP7_FSSA_REFERENCE"] == "auto"


def test_a_restart_keeps_the_reference_it_was_stepped_with():
    from icepack2_tools.fssa import restart_reference_error
    # stepped without the stabilization: either reference may start
    assert restart_reference_error({}, "step", "c.h5") is None
    assert restart_reference_error({}, "start", "c.h5") is None
    # written before the choice existed, so stepped under `start`
    assert restart_reference_error({"fssa_tau": 0.05}, "start", "c.h5") is None
    msg = restart_reference_error({"fssa_tau": 0.05}, "step", "c.h5")
    assert "ISMIP7_FSSA_REFERENCE=start" in msg and "c.h5" in msg
    rec = {"fssa_tau": 0.05, "fssa_reference": "step"}
    assert restart_reference_error(rec, "step", "c.h5") is None
    assert "ISMIP7_FSSA_REFERENCE=step" in restart_reference_error(rec, "start", "c.h5")
