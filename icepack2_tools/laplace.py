r"""Uncertainty from the inferred friction and fluidity.

The inversion returns one state, the MAP, and the Laplace approximation turns
it into a Gaussian posterior: the prior covariance ``Gamma`` updated by the
Gauss-Newton Hessian ``H`` of the data term at the MAP,
``Gamma_post = (H + Gamma^-1)^-1``. ``antarctica/scripts/run_eigendec.py``
computes the leading generalised eigenpairs ``H v_i = lambda_i Gamma^-1 v_i``
(``v_i`` normalised to ``v_i . Gamma^-1 v_j = delta_ij``); in that basis

    Gamma_post = Gamma - sum_i lambda_i / (1 + lambda_i) v_i v_i^T,

so a quantity ``Q`` the model computes from the controls, linearised at the
MAP with gradient ``g = dQ/dm``, has posterior variance

    sigma_post^2 = g . Gamma g - sum_i lambda_i / (1 + lambda_i) (g . v_i)^2.

A truncated set of modes leaves out only reductions, so the variance it
gives is an upper bound on the full Laplace variance.

When the data dominate the prior in more directions than an eigensolve can
hold, ``posterior_variance_cg`` brackets the same variance instead: it solves
``(H + Gamma^-1) x = g`` by conjugate gradients preconditioned with
``Gamma``, where ``g . x_k`` rises to the variance from below and Gauss-Radau
quadrature on the same recurrence bounds it from above. The eigenpairs seed
and deflate that solve, so its zero-iteration upper bound is the truncated
variance above and every Hessian action tightens both sides.

The quantities here are snapshot quantities of the momentum balance at the
MAP state: the ice discharged from grounded ice, across the grounding line
into floating cells and across the mesh boundary, per drainage basin. They
need one adjoint solve each and no time stepping.
"""
import numpy as np


def posterior_variance(prior_variance, projections, eigenvalues):
    r"""Laplace posterior variance of a linearised quantity.

    ``prior_variance`` is ``g . Gamma g``; ``projections`` are ``g . v_i``
    on the ``Gamma^-1``-orthonormal modes and ``eigenvalues`` their
    ``lambda_i``. Returns ``g . Gamma g - sum_i lambda_i/(1+lambda_i)
    (g . v_i)^2``, the full posterior variance when every mode is given and
    an upper bound on it otherwise."""
    p = np.asarray(projections, dtype=float)
    lam = np.asarray(eigenvalues, dtype=float)
    if p.shape != lam.shape:
        raise ValueError(f"{p.size} projections for {lam.size} eigenvalues")
    if np.any(lam < 0.0):
        raise ValueError("a Gauss-Newton Hessian has no negative eigenvalues; "
                         f"got {lam[lam < 0.0][:3]}")
    return float(prior_variance) - float(np.sum(lam / (1.0 + lam) * p ** 2))


def arpack_ncv(k):
    r"""The number of Lanczos vectors ARPACK keeps per restart cycle when none is
    given: scipy's default of ``max(2k + 1, 20)``."""
    return max(2 * k + 1, 20)


def arpack_action_bound(k, ncv, maxiter):
    r"""The most operator actions ARPACK takes for ``k`` pairs with ``ncv``
    Lanczos vectors and ``maxiter`` restarts: one on the start vector, ``ncv``
    for the first factorisation and ``ncv - k`` more per restart."""
    return 1 + ncv + maxiter * (ncv - k)


def leading_eigenpairs(H_action, Ginv_action, G_action, N, k, tol=1e-6, maxiter=None,
                       ncv=None):
    r"""The ``k`` leading generalised eigenpairs ``H v = lambda Gamma^-1 v`` by
    ARPACK, from three actions on global vectors of length ``N``: the Hessian
    ``H``, the prior precision ``Gamma^-1`` and the prior covariance ``Gamma``.

    ``tol`` is ARPACK's relative accuracy on the eigenvalues (0 means machine
    precision). ``maxiter`` bounds the restarts: the solve then takes at most
    ``arpack_action_bound(k, ncv, maxiter)`` Hessian actions, so a job can be
    sized to its wall-clock limit. When ARPACK reaches that bound first, the
    pairs it has converged come back instead of nothing, in which case the
    third value is ``False``. Returns ``(eigenvalues, eigenvectors,
    converged)`` with the eigenvalues descending and the eigenvectors, one per
    column, normalised to ``v_i . Gamma^-1 v_j = delta_ij``; it raises when
    not one pair converged."""
    from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh

    if ncv is None:
        ncv = arpack_ncv(k)
    H_op = LinearOperator((N, N), matvec=H_action, dtype=float)
    Ginv_op = LinearOperator((N, N), matvec=Ginv_action, dtype=float)
    G_op = LinearOperator((N, N), matvec=G_action, dtype=float)
    try:
        vals, vecs = eigsh(H_op, k=k, M=Ginv_op, Minv=G_op, which="LA",
                           tol=tol, ncv=ncv, maxiter=maxiter)
        converged = True
    except ArpackNoConvergence as e:
        vals, vecs = np.asarray(e.eigenvalues, dtype=float), np.asarray(e.eigenvectors, dtype=float)
        converged = False
        if vals.size == 0:
            raise RuntimeError(f"ARPACK converged no eigenpair within {maxiter} restarts "
                               f"({ncv} Lanczos vectors): {e}") from e
    order = np.argsort(-vals)
    return vals[order], vecs[:, order], converged


def cg_jacobi_matrix(alphas, betas):
    r"""The Lanczos matrix ``k`` conjugate-gradient steps define: tridiagonal,
    with ``1/alpha_j + beta_{j-1}/alpha_{j-1}`` on the diagonal and
    ``sqrt(beta_j)/alpha_j`` beside it. Its eigenvalues are the Ritz values of
    the preconditioned operator on the Krylov space those steps built."""
    alphas = np.asarray(alphas, dtype=float)
    betas = np.asarray(betas, dtype=float)
    k = alphas.size
    if betas.size < max(k - 1, 0):
        raise ValueError(f"{k} steps need {k - 1} beta coefficients, got {betas.size}")
    T = np.zeros((k, k))
    for j in range(k):
        T[j, j] = 1.0 / alphas[j] + (betas[j - 1] / alphas[j - 1] if j > 0 else 0.0)
        if j + 1 < k:
            T[j, j + 1] = T[j + 1, j] = np.sqrt(betas[j]) / alphas[j]
    return T


def _radau_upper(alphas, betas, rho0, a):
    r"""Gauss-Radau estimate of ``r_0 . A^-1 r_0`` after ``m`` preconditioned
    conjugate-gradient steps, with one quadrature node fixed at ``a``: the
    ``(m+1) x (m+1)`` Lanczos matrix whose last diagonal entry is moved so
    that ``a`` is an eigenvalue (Golub and Meurant, Matrices, Moments and
    Quadrature, Sect. 6.2), inverted at its first entry. Returns ``nan`` when
    the modified matrix is not positive definite (a Ritz value at or below
    ``a`` through roundoff)."""
    import scipy.linalg

    m = len(alphas)
    d = np.array([1.0 / alphas[j] + (betas[j - 1] / alphas[j - 1] if j > 0 else 0.0)
                  for j in range(m)])
    e = np.array([np.sqrt(betas[j]) / alphas[j] for j in range(m)])
    band = np.zeros((3, m))
    band[1] = d - a
    band[0, 1:] = e[:m - 1]
    band[2, :-1] = e[:m - 1]
    rhs = np.zeros(m)
    rhs[-1] = e[m - 1] ** 2
    try:
        delta = scipy.linalg.solve_banded((1, 1), band, rhs)
    except (ValueError, np.linalg.LinAlgError):
        return np.nan
    omega = a + delta[-1]
    if not np.isfinite(omega) or omega <= 0.0:
        return np.nan
    band = np.zeros((3, m + 1))
    band[1] = np.append(d, omega)
    band[0, 1:] = e
    band[2, :-1] = e
    unit = np.zeros(m + 1)
    unit[0] = 1.0
    try:
        y = scipy.linalg.solve_banded((1, 1), band, unit)
    except (ValueError, np.linalg.LinAlgError):
        return np.nan
    value = rho0 * y[0]
    return value if np.isfinite(value) and value >= 0.0 else np.nan


def quadrature_bounds(alphas, betas, rho0, lower_eigenvalue=1.0):
    r"""Bounds on ``r_0 . A^-1 r_0`` from the coefficients of preconditioned
    conjugate gradients on ``A x = r_0``: ``alphas`` the step lengths,
    ``betas`` the direction updates, ``rho0 = r_0 . M r_0`` with ``M`` the
    preconditioner. Entry ``m`` of each array is the bound after ``m`` steps.

    The lower bounds are the Gauss rule, ``sum_{j<m} alpha_j rho_j``, which is
    the ``A``-norm of the iterate and rises monotonically to the exact value.
    The upper bounds are the Gauss-Radau rule with a node fixed at
    ``lower_eigenvalue``, which must lie at or below the spectrum of the
    preconditioned operator ``M^1/2 A M^1/2``; when roundoff defeats that rule
    at a step the previous upper bound is kept."""
    alphas = np.asarray(alphas, dtype=float)
    betas = np.asarray(betas, dtype=float)
    k = alphas.size
    if betas.size < k:
        raise ValueError(f"{k} steps need {k} beta coefficients, got {betas.size}")
    if lower_eigenvalue <= 0.0:
        raise ValueError(f"the lower eigenvalue must be positive, got {lower_eigenvalue}")
    rhos = rho0 * np.concatenate([[1.0], np.cumprod(betas[:k])])
    lower = np.concatenate([[0.0], np.cumsum(alphas * rhos[:k])])
    upper = np.empty(k + 1)
    upper[0] = rho0 / lower_eigenvalue
    for m in range(1, k + 1):
        value = _radau_upper(alphas[:m], betas[:m], rho0, lower_eigenvalue)
        upper[m] = min(upper[m - 1], value) if np.isfinite(value) else upper[m - 1]
    return lower, upper


def posterior_variance_cg(H_action, Ginv_action, G_action, g, iterations, modes=None,
                          eigenvalues=None, lower_eigenvalue=1.0, state=None, callback=None,
                          stop=None):
    r"""Bracket the Laplace posterior variance ``g . (H + Gamma^-1)^-1 g`` of a
    quantity with gradient ``g`` (a dual vector) by conjugate gradients on
    ``(H + Gamma^-1) x = g`` preconditioned with the prior covariance, from
    the three actions ``leading_eigenpairs`` uses.

    After ``k`` steps ``g . x_k`` is the ``A``-norm of the iterate: it rises
    monotonically to the variance, a lower bound. Gauss-Radau quadrature on
    the same recurrence (``quadrature_bounds``) gives an upper bound, with
    ``lower_eigenvalue`` at or below the spectrum of ``I + Gamma^1/2 H
    Gamma^1/2``: 1 for a positive semi-definite ``H``.

    ``modes`` (one ``Gamma^-1``-orthonormal eigenvector of ``H v = lambda
    Gamma^-1 v`` per row) with their ``eigenvalues`` seed the solve with the
    truncated posterior ``sum_i (g . v_i)^2 / (1 + lambda_i)`` and deflate
    those directions from the iteration. The upper bound at zero iterations
    is then the truncated-spectrum variance of ``posterior_variance``, and
    each iteration tightens both bounds. The seed costs one Hessian action and
    one prior-precision action per mode; each iteration one Hessian action,
    one precision action and one covariance action.

    ``state`` is a returned dictionary to continue from (the same modes must
    be given); ``stop(k)`` ends the loop early when it returns true before
    step ``k``; ``callback(record)`` sees the bounds after each step, and
    after step 0 on a fresh start, with the step histories ``lower_by_step``
    and ``upper_by_step`` and the vectors and coefficients to continue from
    (the keys of ``state``), so a record is itself a state. Returns a
    dictionary with ``lower`` and ``upper`` (entry ``k`` after ``k`` steps,
    including the seed), ``rho``, ``alphas``, ``betas``, ``seed``,
    ``projections``, ``actions`` (Hessian actions taken in this call) and the
    vectors ``x``, ``r``, ``p`` to continue from."""
    g = np.asarray(g, dtype=float)
    a = float(lower_eigenvalue)
    actions = 0
    V = W = None
    if modes is not None and len(modes) > 0:
        V = np.asarray(modes, dtype=float)
        lam = np.asarray(eigenvalues, dtype=float)
        if V.shape != (lam.size, g.size):
            raise ValueError(f"{V.shape[0]} modes of length {V.shape[1]} for "
                             f"{lam.size} eigenvalues and a gradient of length {g.size}")
        W = np.array([np.asarray(Ginv_action(v), dtype=float) for v in V])

    def deflate(r):
        return r - W.T @ (V @ r) if V is not None else r

    def report(step):
        if callback is None:
            return
        lower, upper = quadrature_bounds(alphas, betas, rhos[0], a)
        callback({"iteration": step, "lower": seed + lower[-1], "upper": seed + upper[-1],
                  "lower_by_step": seed + lower, "upper_by_step": seed + upper,
                  "actions": actions, "x": x, "r": r, "p": p, "alphas": list(alphas),
                  "betas": list(betas), "rho": list(rhos), "seed": seed,
                  "projections": projections})

    if state is None:
        if V is not None:
            projections = V @ g
            x = V.T @ (projections / (1.0 + lam))
            Hx = np.asarray(H_action(x), dtype=float)
            actions += 1
            Ax = np.asarray(Ginv_action(x), dtype=float) + Hx
            seed = float(x @ Ax)
            r = deflate(g - Ax)
        else:
            projections = np.zeros(0)
            x = np.zeros_like(g)
            seed = 0.0
            r = g.copy()
        p = np.asarray(G_action(r), dtype=float)
        rho = float(r @ p)
        alphas, betas, rhos = [], [], [rho]
        report(0)
    else:
        x, r, p = (np.array(state[k_], dtype=float) for k_ in ("x", "r", "p"))
        alphas, betas, rhos = list(state["alphas"]), list(state["betas"]), list(state["rho"])
        seed, projections = float(state["seed"]), np.asarray(state["projections"], dtype=float)
        rho = rhos[-1]
    done = len(alphas)
    for k in range(done, done + iterations):
        if (stop is not None and stop(k)) or rho <= 0.0:
            break
        q = np.asarray(Ginv_action(p), dtype=float) + np.asarray(H_action(p), dtype=float)
        actions += 1
        pq = float(p @ q)
        if pq <= 0.0:
            raise RuntimeError(f"conjugate gradients lost positivity at step {k}: "
                               f"p . A p = {pq:.3e}")
        alpha = rho / pq
        x = x + alpha * p
        r = deflate(r - alpha * q)
        z = np.asarray(G_action(r), dtype=float)
        rho_new = float(r @ z)
        beta = rho_new / rho if rho > 0.0 else 0.0
        p = z + beta * p
        alphas.append(alpha)
        betas.append(beta)
        rhos.append(rho_new)
        rho = rho_new
        report(k + 1)
    lower, upper = quadrature_bounds(alphas, betas, rhos[0], a)
    return {"lower": seed + lower, "upper": seed + upper, "rho": rhos, "alphas": alphas,
            "betas": betas, "seed": seed, "projections": projections, "actions": actions,
            "iterations": len(alphas), "x": x, "r": r, "p": p}


def load_modes(path, space, count=None):
    r"""The eigenmodes ``run_eigendec.py`` saved in ``path`` (``mode_<i>_alpha``,
    ``mode_<i>_phi`` on the mesh it saved alongside), as pairs of Functions on
    ``space``, the first ``count`` of them (all when ``None``).

    The saved mesh must be the live mesh (``transfer.meshes_match``). When
    every rank holds its coordinates in the saved order the dof arrays are
    copied; otherwise each mode is interpolated across (same geometry, other
    numbering). Returns ``(modes, copied)``. Collective."""
    import firedrake as fd
    from mpi4py import MPI
    from icepack2_tools.transfer import interpolate_with_fill, meshes_match

    mesh = space.mesh()
    comm = mesh.comm
    with fd.CheckpointFile(path, "r") as chk:
        n_modes = int(chk.get_attr("/", "n_modes"))
        k = n_modes if count is None else min(int(count), n_modes)
        saved = chk.load_mesh()
        if not meshes_match(saved, mesh, comm):
            raise ValueError(f"{path} holds modes of another mesh")
        xa, xb = saved.coordinates.dat.data_ro, mesh.coordinates.dat.data_ro
        same_order = xa.shape == xb.shape and bool(np.allclose(xa, xb))
        copied = comm.allreduce(int(same_order), op=MPI.MIN) == 1
        modes = []
        for i in range(k):
            pair = []
            for control in ("alpha", "phi"):
                f = chk.load_function(saved, f"mode_{i:04d}_{control}")
                target = fd.Function(space, name=f.name())
                if copied:
                    target.dat.data[:] = f.dat.data_ro
                else:
                    interpolate_with_fill(target, f, 0.0, comm)
                pair.append(target)
            modes.append(tuple(pair))
    return modes, copied


def is_dg0(space):
    r"""Whether ``space`` is the DG0 space the finite-volume transport moves
    its thickness in: the discharge forms read cell values and their upwind
    facet fluxes, so they need that space. The check reads the element because
    checkouts name the geometry space differently in a run's context."""
    element = space.ufl_element()
    return element.degree() == 0 and element.family() in ("Discontinuous Lagrange", "DQ", "DG")


def grounded_indicator(Q0, s, b, rho_ratio):
    r"""The DG0 indicator of the cells the discharge counts as grounded.

    A cell is grounded when its surface ``s`` lies above the flotation
    surface of its bed ``b`` (``rho_ratio`` is ``rho_I / rho_W``), or when
    its bed is at or above sea level. The second rule is the ISMIP7 writer's
    (``icepack2_tools.ismip7_output.AnnualOutput``, issue #136): an ice-free
    land cell sits exactly at flotation, so the first rule alone reads it as
    afloat, and ice a land-terminating margin pushes onto bare land would
    count as grounding-line flux. With it, the grounding line separates
    grounded cells from marine cells that are not grounded (floating ice or
    open water), and flow onto land books nowhere."""
    import firedrake as fd

    s_float = fd.Function(Q0).interpolate(
        b + fd.max_value(-b, fd.Constant(0.0)) / rho_ratio)
    bed = fd.Function(Q0).interpolate(b)
    grounded = fd.Function(Q0, name="grounded")
    grounded.dat.data[:] = ((s.dat.data_ro - s_float.dat.data_ro > 0.0)
                            | (bed.dat.data_ro >= 0.0)).astype(float)
    return grounded


def grounded_discharge_forms(u, h, grounded, weight=None):
    r"""The ice leaving grounded cells, as two UFL forms in m3/yr of ice.

    ``u`` is the velocity, ``h`` the DG0 thickness, ``grounded`` a DG0
    indicator of the grounded cells (``grounded_indicator``, which counts
    every cell on a bed at or above sea level as grounded) and ``weight`` an
    optional DG0 weight on the grounded side (a basin indicator), read where
    the ice leaves.

    The first form is the grounding-line flux: across every facet between a
    grounded cell and a marine cell that is not grounded, the upwind flux the
    DG0 transport moves, signed positive from grounded to floating, as the
    ISMIP7 writer books ``ligroundf``. Ice flowing back onto a pinning point counts negative, so
    an ice rumple's throughput cancels. The second is the upwind outflow of
    grounded cells across the mesh boundary, which a mesh whose edge is the
    ice front carries at grounded margins."""
    import firedrake as fd

    mesh = h.function_space().mesh()
    n = fd.FacetNormal(mesh)
    un = fd.dot(u, n)
    out = (un + abs(un)) / 2
    flux = out("+") * h("+") - out("-") * h("-")       # upwind, "+" to "-"
    w = weight if weight is not None else fd.Constant(1.0)
    w_plus = w("+") if weight is not None else w
    w_minus = w("-") if weight is not None else w
    gl = (flux * grounded("+") * (1 - grounded("-")) * w_plus
          - flux * grounded("-") * (1 - grounded("+")) * w_minus) * fd.dS
    edge = out * h * grounded * w * fd.ds
    return gl, edge
