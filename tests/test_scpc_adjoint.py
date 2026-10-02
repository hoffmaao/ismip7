r"""The inversion's taped solve under each ISMIP7_INVERSION_LINEAR_SOLVER
mode gives one objective and one gradient.

A 30 km x 10 km grounded slab with a floating tongue (Budd law, n = 3), the
real CG1 x DG0-symmetric-tensor x DG0-vector layout and the inversion's two
residual builders, no data files: the cell-wise law (``build_rc_residual``,
no friction on the tongue) and the sub-element grounding scheme
(``build_subelement_residual``, with the exact front push on an ice-free strip
past 27 km and the friction ``C_ref exp(theta)`` of the exp control). The
gradient of a velocity misfit with respect to both controls is computed
through ``taped_state_solve`` and tlm_adjoint:

- ``scpc_mumps`` condenses exactly, so its gradient is the assembled
  ``full_mumps`` gradient to roundoff (3e-14 to 3e-13 measured);
- ``scpc_gamg`` stops its Krylov solves at the inversion's relative
  tolerance (1e-8 or better measured);
- a Taylor test of the ``scpc_gamg`` gradient converges at second order;
- ``scpc_gamg`` takes its Newton steps whole: 4 iterations here, as the exact
  condensation does. Backtracking on ||F|| (``bt``) took 113 and 20 on these
  slabs, with an exact LU as with GAMG.

Before ``with_quadrature_degree`` both scpc gradients sat 1e-3 off the
assembled one: tlm_adjoint's matrix-free adjoint solve passes no form compiler
parameters, so it differentiated a degree-26 discretization of a degree-4
forward.
"""
import pytest

firedrake = pytest.importorskip("firedrake")

# dual_friction next: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
from icepack2_tools.dual_friction import build_rc_residual     # noqa: E402

from firedrake import (                                         # noqa: E402
    Constant,
    FiniteElement,
    Function,
    FunctionSpace,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    RectangleMesh,
    SpatialCoordinate,
    TensorFunctionSpace,
    VectorFunctionSpace,
    as_vector,
    conditional,
    derivative,
    dx,
    exp,
    inner,
    max_value,
    split,
)
from tlm_adjoint.firedrake import (                             # noqa: E402
    Functional,
    compute_gradient,
    reset_manager,
    start_manager,
    stop_manager,
    taylor_test,
)
from ufl import adjoint                                         # noqa: E402

from icepack2_tools.preconditioners import with_scpc_blocks     # noqa: E402
from icepack2_tools.subelement import (                         # noqa: E402
    build_subelement_residual,
    ice_indicator,
    subelement_from_geometry,
)
from icepack2_tools.solverconfig import (                       # noqa: E402
    inversion_adjoint_parameters,
    inversion_state_parameters,
)
from icepack2_tools.taped_solve import (                        # noqa: E402
    taped_state_solve,
    with_quadrature_degree,
)

RHO_I, RHO_W = 917.0, 1024.0
FCP = {"quadrature_degree": 4}
MODES = ("full_mumps", "scpc_mumps", "scpc_gamg")


@pytest.fixture(scope="module", params=["cellwise", "subelement"])
def slab(request):
    stop_manager()
    law = request.param
    mesh = RectangleMesh(30, 10, 30e3, 10e3)
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Z = V * TensorFunctionSpace(mesh, dg0, symmetry=True) * VectorFunctionSpace(mesh, dg0)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, y = SpatialCoordinate(mesh)
    if law == "cellwise":
        # grounded for x < 20 km, floating beyond
        H = Function(Q0).interpolate(conditional(x < 20e3, 1500.0 - 0.03 * x, 400.0))
        b = Function(Q0).interpolate(conditional(x < 20e3, -200.0 - 0.01 * x, -1500.0))
        s = Function(Q0).interpolate(max_value(b + H, (1.0 - RHO_I / RHO_W) * H))
    else:
        # a thinning wedge on a deepening bed: the grounding line crosses
        # cells, and the ice ends at a front inside the mesh at 27 km
        H = Function(Q0).interpolate(conditional(x < 27e3, 1200.0 - 0.03 * x, 0.0))
        b = Function(Q0).interpolate(-100.0 - 0.04 * x)
        haf = H - Constant(RHO_W / RHO_I) * conditional(b < 0.0, -b, 0.0)
        s = Function(Q0).interpolate(
            conditional(haf > 0.0, b + H, (1.0 - RHO_I / RHO_W) * H))
    f = dict(
        law=law, mesh=mesh, x=x, y=y, Q=Q, Z=Z, H=H, b=b, s=s,
        C_w0=Function(Q).interpolate(
            0.05 + 0.03 * (1.0 + y / 10e3) * (1.0 - x / 40e3)),
        A_prior=Function(Q).interpolate(Constant(20.0) + 5.0 * x / 30e3),
        n_flow=Constant(3.0),
        u_obs=Function(V).interpolate(as_vector([100.0 + 0.02 * x, 2.0 * y / 10e3])),
        theta0=Function(Q, name="theta").interpolate(0.2 * x / 30e3 - 0.1),
        phi0=Function(Q, name="phi").interpolate(0.1 * y / 10e3),
    )
    if law == "subelement":
        f["subelement"] = subelement_from_geometry(mesh, H, b, ice=ice_indicator(H, 1.0))
    # A converged state at the reference controls: climb n from 1, untaped.
    z0 = Function(Z)
    z0.subfunctions[0].interpolate(0.1 * f["u_obs"])
    params = inversion_state_parameters("full_mumps")
    for n in (1.0, 1.5, 2.0, 2.5, 3.0):
        f["n_flow"].assign(n)
        NonlinearVariationalSolver(
            NonlinearVariationalProblem(
                _residual(f, z0, f["theta0"], f["phi0"], scpc=False), z0,
                form_compiler_parameters=FCP),
            solver_parameters=params,
        ).solve()
    f["z0"] = z0
    return f


def _residual(f, z, theta, phi, *, scpc):
    if f["law"] == "cellwise":
        F = build_rc_residual(
            z, theta, phi, H=f["H"], s=f["s"], b=f["b"], C_w0=f["C_w0"],
            A4_base=f["A_prior"], n_flow=f["n_flow"], n_flow_val=3.0,
            m_slide=3.0, tau_c=Constant(0.1), alpha=Constant(1e-2),
            H_ref=Constant(100.0), fric_law="budd", N_ref=None,
            nhat_floor=0.02, nhat_cap=3.0, alpha_gl=0.5, c_w0_floor=0.0,
            h_visc_floor=10.0, ocean_drag=0.0, k_lim=0.0,
        )
    else:
        # the exp control as the inversion passes it: a zero log deviation
        # and the friction C_ref exp(theta) outright
        F = build_subelement_residual(
            z, Constant(0.0), phi, H=f["H"], s=f["s"], b=f["b"],
            C_w0=Constant(0.05) * exp(theta), A4_base=f["A_prior"],
            n_flow=f["n_flow"], n_flow_val=3.0, m_slide=3.0,
            tau_c=Constant(0.1), alpha=Constant(1e-2), H_ref=Constant(100.0),
            subelement=f["subelement"], fric_law="budd", nhat_cap=3.0,
            alpha_gl=0.5, c_w0_floor=0.0, h_visc_floor=10.0, ocean_drag=1e-2,
            h_ocean=10.0, u_lim=0.0, k_lim=0.0, exact_front=True,
        )
    return with_scpc_blocks(F, z) if scpc else F


def _misfit(f, mode, theta, phi):
    z = Function(f["Z"]).assign(f["z0"])
    params = inversion_state_parameters(mode)
    work = taped_state_solve(
        _residual(f, z, theta, phi, scpc=mode != "full_mumps"), z, mode,
        params, inversion_adjoint_parameters(params),
        form_compiler_parameters=FCP,
    )
    u = split(z)[0]
    J = Functional(name="J")
    J.assign(0.5 * inner(u - f["u_obs"], u - f["u_obs"]) * dx)
    return J, work


def _controls(f):
    # off the reference, so the taped solve does real Newton work from z0
    theta = Function(f["Q"], name="theta").interpolate(f["theta0"] + 0.05)
    phi = Function(f["Q"], name="phi").assign(f["phi0"])
    return theta, phi


@pytest.fixture(scope="module")
def gradients(slab):
    out = {}
    for mode in MODES:
        theta, phi = _controls(slab)
        reset_manager()
        start_manager()
        J, work = _misfit(slab, mode, theta, phi)
        stop_manager()
        dJ = compute_gradient(J, [theta, phi])
        out[mode] = dict(J=float(J.value), dJ=dJ, work=work, controls=(theta, phi))
    reset_manager()
    return out


def _rel(a, b):
    with a.dat.vec_ro as va, b.dat.vec_ro as vb:
        d = va.copy()
        d.axpy(-1.0, vb)
        return d.norm() / vb.norm()


@pytest.mark.parametrize("mode, tol", [("scpc_mumps", 1e-10), ("scpc_gamg", 1e-6)])
def test_scpc_gradient_is_the_assembled_gradient(gradients, mode, tol):
    ref, got = gradients["full_mumps"], gradients[mode]
    assert got["J"] == pytest.approx(ref["J"], rel=1e-10)
    for g, g_ref in zip(got["dJ"], ref["dJ"]):
        assert _rel(g, g_ref) < tol


@pytest.mark.parametrize("mode", ["scpc_mumps", "scpc_gamg"])
def test_the_recorded_solve_only_confirms(gradients, mode):
    work = gradients[mode]["work"]
    assert work["converged_reason"] > 0
    assert work["snes_iterations"] > 0
    assert work["condensed_solves"] >= work["linear_iterations"] > 0
    # the tape holds the state the untaped solve reached
    assert work["confirm_step"] == 0.0


def test_gamg_takes_newton_steps_whole(gradients):
    # the exact condensation under NLEQ-ERR takes 4 on both slabs; bt took
    # 113 and 20, under either solver
    assert gradients["scpc_gamg"]["work"]["snes_iterations"] <= 8


def test_scpc_gamg_gradient_passes_a_taylor_test(slab, gradients):
    g = gradients["scpc_gamg"]
    theta, phi = g["controls"]

    def forward(th):
        J, _ = _misfit(slab, "scpc_gamg", th, phi)
        return J

    dM = Function(slab["Q"]).interpolate(0.05 * slab["x"] / 30e3)
    order = taylor_test(forward, theta, J_val=g["J"], dJ=g["dJ"][0], dM=dM,
                        seed=1e-1, size=4)
    assert order > 1.9


def test_the_degree_survives_into_the_adjoint_form(slab):
    f = slab
    z = Function(f["Z"]).assign(f["z0"])
    F = with_quadrature_degree(
        _residual(f, z, f["theta0"], f["phi0"], scpc=True), FCP)
    for form in (F, derivative(F, z), adjoint(derivative(F, z))):
        degrees = {i.metadata().get("quadrature_degree") for i in form.integrals()}
        assert degrees == {4}
