r"""The inversion's taped momentum solve, under each inversion solver mode
(``solverconfig.inversion_solver_mode``, ``ISMIP7_INVERSION_LINEAR_SOLVER``).

tlm_adjoint records the solve as ``EquationSolver(F == 0, z)`` and computes
the adjoint by assembling ``adjoint(derivative(F, z))`` at the recorded state
and solving that form with the adjoint options. It never applies a transpose,
so the SCPC preconditioner, which has none, serves the adjoint as it serves
the forward.

How the forward gets to that state differs. ``full_mumps`` lets the recorded
equation solve itself, as it always has. Under ``scpc_*`` the Newton solve
needs a Jacobian frozen at the iterate SCPC assembled its condensed system at
(``preconditioners.frozen_linearization``), and the frozen Jacobian is
refreshed by a ``pre_jacobian_callback`` that tlm_adjoint's patched
``NonlinearVariationalSolver.solve`` refuses while it annotates ("Callbacks
not supported"). The Newton solve therefore runs with the manager paused, and
the recorded equation then confirms the converged state: its absolute
tolerance is set just above the residual the paused solve reached, so SNES
accepts the state at iteration 0, with no linear solve, and the tape holds
``F`` and the converged ``z`` exactly as a recorded solve would. The driver
already writes ``z`` outside the tape before a recorded solve the same way
(the startup ramp, and the re-ramp rescue at a failed trial point).
"""

import ufl
from firedrake import (
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    assemble,
)
from tlm_adjoint.firedrake import EquationSolver, paused_manager

from icepack2_tools.preconditioners import frozen_linearization

# The recorded confirmation re-evaluates the residual at the state the paused
# solve left, which on one partition is the same number; the margin only
# keeps a last-digit difference from turning the confirmation into a solve.
CONFIRM_ATOL_MARGIN = 1.01


def with_quadrature_degree(F, form_compiler_parameters):
    r"""``F`` with the quadrature degree of ``form_compiler_parameters``
    written into each integral's metadata, where every form derived from it
    (its Jacobian, and the adjoint of that) inherits it.

    Under ``mat_type: matfree`` tlm_adjoint solves the adjoint with
    ``solve(adjoint(J) == b)`` and no form compiler parameters, so that
    operator and the SCPC condensation of it are assembled at UFL's estimated
    degree (up to 26 here) while the forward ran at 4: the gradient of another
    discretization, 1e-3 off the assembled path's on a synthetic slab
    (tests/test_scpc_adjoint.py). An integral that names its own degree keeps
    it, as it would over the parameters."""
    degree = (form_compiler_parameters or {}).get("quadrature_degree")
    if degree is None:
        return F
    if not isinstance(F, ufl.Form):
        raise TypeError(f"expected a ufl.Form, not {type(F).__name__}")
    return ufl.Form([
        integral.reconstruct(
            metadata={"quadrature_degree": degree, **integral.metadata()})
        for integral in F.integrals()
    ])


def _condensed_work(snes):
    r"""Solves and iterations on the condensed velocity system, from the
    ISMIP7SCPC counters (SNES's own linear count misses the line search's)."""
    try:
        ctx = snes.getKSP().getPC().getPythonContext()
    except Exception:
        return None, None
    return (getattr(ctx, "condensed_solves", None),
            getattr(ctx, "condensed_iterations", None))


def taped_state_solve(F, z, mode, params, adjoint_params, *,
                      form_compiler_parameters=None,
                      options_prefix="ismip7_inversion_state_"):
    r"""Solve ``F(z) = 0`` and record it on the tlm_adjoint tape.

    ``params`` are the forward's options under ``mode``
    (``solverconfig.inversion_state_parameters``) and ``adjoint_params`` the
    adjoint's (``solverconfig.inversion_adjoint_parameters``). Under
    ``scpc_*``, ``F`` must carry the SCPC structural-zero blocks
    (``preconditioners.with_scpc_blocks``). A solve that does not converge
    raises ``firedrake.ConvergenceError``, as the recorded solve does.

    Returns a dict of the untaped solve's work under ``scpc_*`` (empty under
    ``full_mumps``, whose work happens inside the recorded solve).
    """
    if mode == "full_mumps":
        EquationSolver(
            F == 0,
            z,
            solver_parameters=params,
            adjoint_solver_parameters=adjoint_params,
            form_compiler_parameters=form_compiler_parameters,
        ).solve()
        return {}

    F = with_quadrature_degree(F, form_compiler_parameters)
    J, pre_jacobian = frozen_linearization(F, z)
    with paused_manager():
        solver = NonlinearVariationalSolver(
            NonlinearVariationalProblem(
                F, z, J=J, form_compiler_parameters=form_compiler_parameters
            ),
            solver_parameters=params,
            options_prefix=options_prefix,
            pre_jacobian_callback=pre_jacobian,
        )
        solver.solve()
        snes = solver.snes
        condensed_solves, condensed_iterations = _condensed_work(snes)
        work = {
            "snes_iterations": int(snes.getIterationNumber()),
            "linear_iterations": int(snes.getLinearSolveIterations()),
            "converged_reason": int(snes.getConvergedReason()),
            "condensed_solves": condensed_solves,
            "condensed_iterations": condensed_iterations,
        }
        with assemble(
            F, form_compiler_parameters=form_compiler_parameters
        ).dat.vec_ro as residual:
            fnorm = float(residual.norm())
    work["fnorm"] = fnorm

    confirm = dict(params)
    confirm["snes_atol"] = max(
        float(params.get("snes_atol", 0.0)), CONFIRM_ATOL_MARGIN * fnorm
    )
    with z.dat.vec_ro as v:
        before = v.copy()
    EquationSolver(
        F == 0,
        z,
        solver_parameters=confirm,
        adjoint_solver_parameters=adjoint_params,
        form_compiler_parameters=form_compiler_parameters,
    ).solve()
    # Zero unless the confirmation took a Newton step: the tape then holds a
    # state the paused solve did not reach, still a converged one.
    with z.dat.vec_ro as v:
        before.axpy(-1.0, v)
    work["confirm_step"] = float(before.norm())
    before.destroy()
    return work
