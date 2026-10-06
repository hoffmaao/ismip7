r"""The inversion's taped momentum solve, under each inversion solver mode
(``solverconfig.inversion_solver_mode``, ``ISMIP7_INVERSION_LINEAR_SOLVER``).

tlm_adjoint records the solve as ``EquationSolver(F == 0, z)`` and computes
the adjoint by assembling ``adjoint(derivative(F, z))`` at the recorded state
and solving that form with the adjoint options. It never applies a transpose,
so the SCPC preconditioner, which has none, serves the adjoint as it serves
the forward.

How the forward gets to that state differs. With the direct forward disabled,
``full_mumps`` lets the recorded equation solve itself, as it always has.
Under ``scpc_*`` the Newton solve needs a Jacobian frozen at the iterate SCPC
assembled its condensed system at (``preconditioners.frozen_linearization``),
and the frozen Jacobian is refreshed by a ``pre_jacobian_callback`` that
tlm_adjoint's patched ``NonlinearVariationalSolver.solve`` refuses while it
annotates ("Callbacks not supported"). The Newton solve therefore runs with
the manager paused, and the recorded equation then confirms the converged
state: its absolute tolerance is set just above the residual the paused solve
reached, so SNES accepts the state at iteration 0, with no linear solve, and
the tape holds ``F`` and the converged ``z`` exactly as a recorded solve
would. The driver already writes ``z`` outside the tape before a recorded
solve the same way (the startup ramp, and the re-ramp rescue at a failed
trial point).

The paused solver's Firedrake context can outlive the call. In the condensed
modes that context includes the frozen linearization's state Function and,
from its first linear solve, the SCPC context that holds the GAMG hierarchy.
``StateSolverCache`` keeps it for the next call with the same form. The
recorded confirmation and the adjoint are tlm_adjoint's own solves and build
theirs each time.

The direct forward (``ISMIP7_DIRECT_FORWARD``, on by default) takes the same
route under every mode, ``full_mumps`` included: one paused Newton solve at
the full exponents from the state ``z`` holds, under
``solverconfig.direct_forward_parameters``, then the recorded confirmation.
A direct solve that does not converge leaves ``z`` at its entry state and
raises, so the driver's failed-trial rescue starts from the last converged
state. ``full_mumps`` keeps its Newton Jacobian live; the condensed modes
freeze it as above.

Every solver this module builds and drops is destroyed as it is dropped
(``release_solver``, issue #161): the recorded equations'
(``ReleasingEquationSolver``: the Newton solver of the recorded solve, whose
assembled Jacobian is allocated even when SNES exits at iteration 0, and under
``full_mumps`` the adjoint's assembled operator, solver and LU) and the paused
solver a new form or a failed solve replaces. Left to themselves they are
freed late, and some never. A Firedrake ``NonlinearVariationalSolver`` sits
in a reference cycle, so it and its Jacobian live until Python's cyclic
collector runs, and on more than one rank petsc4py then hands its PETSc
objects to the next ``PetscGarbageCleanup``. Under ``full_mumps`` on the 2 km
mesh the resident set climbed 8 to 20 MiB a rank an evaluation that way, and
28 to 37 matrices were still alive at exit against 10 with the release
(README "Inversion solver", issue #161 run records).
"""

from time import perf_counter

import ufl
from firedrake import (
    ConvergenceError,
    Function,
    NonlinearVariationalProblem,
    NonlinearVariationalSolver,
    assemble,
)
from firedrake.utils import ScalarType
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
    ISMIP7SCPC counters (SNES's own linear count misses the line search's).
    The counters run from the PC's first setup, so a reused solver's work
    is the difference across the call. ``None`` under any other PC: petsc4py's
    ``getPythonContext`` on a PC of another type (the direct forward's MUMPS
    LU under ``full_mumps``) is a segmentation fault, not an exception."""
    pc = snes.getKSP().getPC()
    if pc.getType() != "python":
        return None, None
    try:
        ctx = pc.getPythonContext()
    except Exception:
        return None, None
    return (getattr(ctx, "condensed_solves", None),
            getattr(ctx, "condensed_iterations", None))


def _since(after, before):
    return None if after is None else after - (before or 0)


def release_solver(solver):
    r"""Destroy what a Firedrake variational solver built, now: its SNES with
    the KSP, PC and any factor, its context's Jacobian and preconditioning
    matrices (a matrix-free one's Python context with them), its work vector,
    and a ``LinearSolver``'s operator. Collective: every rank builds the same
    solvers, so every rank releases them in the same order. The solver cannot
    solve again."""
    ctx = getattr(solver, "_ctx", None)
    # cached properties, present only once the solver has made them
    mats = [ctx.__dict__.get(name) for name in ("_jac", "_pjac")] if ctx else []
    mats.append(getattr(solver, "A", None))
    solver.snes.destroy()
    for mat in mats:
        petscmat = getattr(mat, "petscmat", None)
        if petscmat is not None:
            petscmat.destroy()
    work = getattr(solver, "_work", None)
    if work is not None:
        work.destroy()


class ReleasingEquationSolver(EquationSolver):
    r"""tlm_adjoint's ``EquationSolver``, releasing every solver it builds for
    one solve (``release_solver``) once that solve returns.

    Those are the forward's Newton solver, which tlm_adjoint builds through
    ``firedrake.solve`` and is built here the same way so it can be kept for
    release, and the assembled operator and ``LinearSolver`` of an uncached
    linear forward or adjoint solve. Cached solvers stay tlm_adjoint's, and
    so does the matrix-free adjoint (``ISMIP7SCPC.destroy`` covers its
    condensed solver). Nothing numerical differs from the base class: the
    inversion's objectives and gradients are the same to every printed digit
    (issue #161 run records)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._built = []

    def _assemble_linear_solver(self, *args, **kwargs):
        solver, b_bc = super()._assemble_linear_solver(*args, **kwargs)
        self._built.append(solver)
        return solver, b_bc

    def _release_built(self):
        while self._built:
            release_solver(self._built.pop(0))

    def forward_solve(self, x, deps=None):
        try:
            if (self._linear
                    and self._solver_parameters.get("mat_type", "aij") != "matfree"):
                super().forward_solve(x, deps)
            else:
                self._newton_solve(x, deps)
        finally:
            self._release_built()

    def adjoint_jacobian_solve(self, adj_x, nl_deps, b):
        try:
            return super().adjoint_jacobian_solve(adj_x, nl_deps, b)
        finally:
            self._release_built()

    def _newton_solve(self, x, deps):
        # EquationSolver.forward_solve's solve(F == 0, ...) as
        # tlm_adjoint.firedrake.backend_interface.solve and
        # firedrake.solving._solve_varproblem carry it out, keeping the solver.
        # Like them, this adds scalar_type to the equation's own form
        # compiler parameters, which the adjoint then assembles with.
        bcs = tuple(self._reconstruct_bcs(deps=deps))
        J = self._J if self._nl_solve_J is None else self._nl_solve_J
        params = dict(self._solver_parameters)
        extra = params.pop("tlm_adjoint", {})
        if "pre_apply_bcs" in extra:
            raise TypeError("Cannot pass both pre_apply_bcs argument and "
                            "solver parameter")
        fcp = self._form_compiler_parameters
        fcp["scalar_type"] = ScalarType
        problem = NonlinearVariationalProblem(
            self._replace(self._F, deps), x, bcs, self._replace(J, deps), None,
            form_compiler_parameters=fcp, restrict=False)
        solver = NonlinearVariationalSolver(
            problem, solver_parameters=params,
            nullspace=extra.get("nullspace"),
            transpose_nullspace=extra.get("transpose_nullspace"),
            near_nullspace=extra.get("near_nullspace"),
            options_prefix=extra.get("options_prefix"),
            appctx={}, pre_apply_bcs=True)
        try:
            solver.solve()
        finally:
            release_solver(solver)


class StateSolverCache:
    r"""The paused Newton solver of the last residual form
    ``taped_state_solve`` was given, kept for the next call with that form.

    A call reuses it when its form is the same object, with the same state
    Function, options, form compiler parameters and options prefix; any other
    call builds a new solver and drops this one, so a caller whose form
    changes on every call holds one solver at a time. The inversion builds
    its taped form once for the run (``inversion_icepack2.py``,
    ``_residual_at``), so one solver serves every evaluation.

    A reused solver reaches the state a new one does, to roundoff. In a
    condensed mode each Jacobian writes its iterate into the frozen
    linearization's state before anything reads it; NLEQ-ERR clears its step
    history at iteration 0 of every SNES solve; SCPC reassembles its condensed
    matrix on every new Jacobian and starts its Krylov solves from zero. What
    the condensed preconditioner built at its first setup stays, as it does
    across the Newton iterations of one solve and across every step of a
    transient run, which keeps one solver throughout: GAMG redoes its Galerkin
    products and Chebyshev eigenvalue estimates on its first interpolation
    (``pc_gamg_reuse_interpolation`` and ``pc_gamg_recompute_esteig``, both
    true by default), and an LU refactors on its first symbolic analysis
    (README, "Inversion solver", for what that changed). A solve that raises
    drops the solver, so a retry starts from a new one, as it did before the
    cache. A solver dropped for any reason is released (``release_solver``).
    """

    def __init__(self):
        self._entry = None

    def clear(self):
        r"""Release the kept solver, if any. Collective."""
        entry, self._entry = self._entry, None
        if entry is not None:
            release_solver(entry["solver"])

    def get(self, F, z, params, form_compiler_parameters, options_prefix,
            *, frozen=True):
        r"""``(F_q, solver, reused)``: ``F`` with the quadrature degree in its
        integrals (``with_quadrature_degree``) and the Newton solver of
        ``F_q(z) = 0``, with its Jacobian frozen at the iterate (``frozen``,
        the condensed modes) or live (``full_mumps``). Call it with the
        manager paused."""
        key = (dict(params), dict(form_compiler_parameters or {}),
               options_prefix, bool(frozen))
        entry = self._entry
        if (entry is not None and entry["F"] is F and entry["z"] is z
                and entry["key"] == key):
            return entry["F_q"], entry["solver"], True
        self.clear()
        F_q = with_quadrature_degree(F, form_compiler_parameters)
        J, pre_jacobian = (frozen_linearization(F_q, z) if frozen
                           else (None, None))
        solver = NonlinearVariationalSolver(
            NonlinearVariationalProblem(
                F_q, z, J=J, form_compiler_parameters=form_compiler_parameters
            ),
            solver_parameters=params,
            options_prefix=options_prefix,
            pre_jacobian_callback=pre_jacobian,
        )
        self._entry = {"F": F, "z": z, "key": key, "F_q": F_q,
                       "solver": solver, "residual": None, "z_entry": None}
        return F_q, solver, False

    def save_entry_state(self):
        r"""Copy ``z`` into a buffer kept with the solver, for
        :meth:`restore_entry_state` after a direct solve that fails."""
        entry = self._entry
        if entry["z_entry"] is None:
            entry["z_entry"] = Function(entry["z"].function_space())
        entry["z_entry"].assign(entry["z"])
        return entry["z_entry"]

    def residual_norm(self, form_compiler_parameters):
        r"""``||F_q(z)||`` of the cached form, assembled into one buffer."""
        entry = self._entry
        entry["residual"] = assemble(
            entry["F_q"], tensor=entry["residual"],
            form_compiler_parameters=form_compiler_parameters)
        with entry["residual"].dat.vec_ro as residual:
            return float(residual.norm())


def taped_state_solve(F, z, mode, params, adjoint_params, *,
                      form_compiler_parameters=None,
                      options_prefix="ismip7_inversion_state_",
                      cache=None, direct=False):
    r"""Solve ``F(z) = 0`` and record it on the tlm_adjoint tape.

    ``params`` are the forward's options under ``mode``
    (``solverconfig.inversion_state_parameters``) and ``adjoint_params`` the
    adjoint's (``solverconfig.inversion_adjoint_parameters``). Under
    ``scpc_*``, ``F`` must carry the SCPC structural-zero blocks
    (``preconditioners.with_scpc_blocks``), and ``cache``, a
    ``StateSolverCache``, keeps the paused solver for the next call with the
    same ``F``; without one each call builds its own and releases it before
    returning. A solve that does not converge raises
    ``firedrake.ConvergenceError``, as the recorded solve does.

    With ``direct`` every mode, ``full_mumps`` included, solves untaped under
    ``params``, which are then the direct solve's options
    (``solverconfig.direct_forward_parameters``), and a solve that does not
    converge leaves ``z`` at its entry state before it raises, with its SNES
    reason, iterations and residual in the message.

    Returns a dict of the untaped solve's work, with ``reused`` saying
    whether its solver came from ``cache`` (empty under ``full_mumps``
    without ``direct``, whose work happens inside the recorded solve). Its
    ``seconds`` member is this rank's untaped solve time. The inversion timing
    record separately reduces the surrounding forward time across ranks.
    """
    if mode == "full_mumps" and not direct:
        ReleasingEquationSolver(
            F == 0,
            z,
            solver_parameters=params,
            adjoint_solver_parameters=adjoint_params,
            form_compiler_parameters=form_compiler_parameters,
        ).solve()
        return {}

    own_cache = cache is None
    if own_cache:
        cache = StateSolverCache()
    with paused_manager():
        F, solver, reused = cache.get(
            F, z, params, form_compiler_parameters, options_prefix,
            frozen=mode != "full_mumps")
        snes = solver.snes
        solves_before, iterations_before = _condensed_work(snes)
        z_entry = cache.save_entry_state() if direct else None
        t0 = perf_counter()
        try:
            solver.solve()
        except ConvergenceError as err:
            # read before the release destroys the SNES
            failure = (f"direct forward did not converge: SNES reason "
                       f"{int(snes.getConvergedReason())}, "
                       f"{int(snes.getIterationNumber())} iterations, "
                       f"||F|| {snes.getFunctionNorm():.3e}")
            cache.clear()
            if z_entry is None:
                raise
            z.assign(z_entry)
            raise ConvergenceError(failure) from err
        except BaseException:
            cache.clear()
            raise
        seconds = perf_counter() - t0
        condensed_solves, condensed_iterations = _condensed_work(snes)
        work = {
            "snes_iterations": int(snes.getIterationNumber()),
            "linear_iterations": int(snes.getLinearSolveIterations()),
            "converged_reason": int(snes.getConvergedReason()),
            "condensed_solves": _since(condensed_solves, solves_before),
            "condensed_iterations": _since(
                condensed_iterations, iterations_before),
            "reused": reused,
            "seconds": seconds,
        }
        fnorm = cache.residual_norm(form_compiler_parameters)
        if own_cache:
            cache.clear()
    work["fnorm"] = fnorm

    confirm = dict(params)
    confirm["snes_atol"] = max(
        float(params.get("snes_atol", 0.0)), CONFIRM_ATOL_MARGIN * fnorm
    )
    with z.dat.vec_ro as v:
        before = v.copy()
    ReleasingEquationSolver(
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
