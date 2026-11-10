r"""The Laplace posterior of a quantity computed from the inferred controls
(icepack2_tools.laplace), and the grounded-discharge forms it is applied to.

The variance formula is checked against the exact posterior of a small
linear-Gaussian problem, where ``Gamma_post = (G^T W G + Gamma^-1)^-1``
holds without approximation."""
import numpy as np
import pytest
import scipy.linalg

from icepack2_tools.laplace import posterior_variance


def _linear_gaussian(n=7, m=4, seed=3):
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(n, n))
    prior_cov = A @ A.T + n * np.eye(n)
    G = rng.normal(size=(m, n))
    W = np.diag(rng.uniform(0.5, 4.0, size=m))
    H = G.T @ W @ G
    lam, V = scipy.linalg.eigh(H, np.linalg.inv(prior_cov))   # V^T Gamma^-1 V = I
    g = rng.normal(size=n)
    exact = g @ np.linalg.solve(H + np.linalg.inv(prior_cov), g)
    return prior_cov, lam, V, g, exact


def test_every_mode_gives_the_exact_linear_gaussian_variance():
    prior_cov, lam, V, g, exact = _linear_gaussian()
    lam = np.clip(lam, 0.0, None)            # roundoff on the null space
    got = posterior_variance(g @ prior_cov @ g, V.T @ g, lam)
    assert got == pytest.approx(exact, rel=1e-10)


def test_leading_modes_alone_bound_the_variance_from_above():
    prior_cov, lam, V, g, exact = _linear_gaussian()
    order = np.argsort(-lam)
    previous = g @ prior_cov @ g
    for k in range(1, 4):
        top = order[:k]
        got = posterior_variance(g @ prior_cov @ g, V[:, top].T @ g, np.clip(lam[top], 0.0, None))
        assert exact - 1e-12 <= got <= previous + 1e-12
        previous = got


def test_a_negative_eigenvalue_is_refused():
    with pytest.raises(ValueError, match="negative"):
        posterior_variance(1.0, [0.1, 0.2], [2.0, -1.0])


def test_projections_and_eigenvalues_must_pair_up():
    with pytest.raises(ValueError, match="projections"):
        posterior_variance(1.0, [0.1, 0.2, 0.3], [2.0, 1.0])


def _diagonal_problem(spectrum, seed=5):
    """A generalised eigenproblem H v = lambda Gamma^-1 v with diagonal H and
    Gamma, built to have the given eigenvalues: its eigenvectors are the unit
    vectors scaled to v . Gamma^-1 v = 1, and every action is counted."""
    from icepack2_tools.laplace import leading_eigenpairs
    lam = np.asarray(spectrum, dtype=float)
    rng = np.random.default_rng(seed)
    cov = rng.uniform(0.5, 2.0, size=lam.size)      # Gamma
    hess = lam / cov                                 # Gamma H = diag(lam)
    counts = {"H": 0}

    def H_action(x):
        counts["H"] += 1
        return hess * x

    def solve(k, **kw):
        return leading_eigenpairs(H_action, lambda x: x / cov, lambda x: cov * x, lam.size, k, **kw)

    return solve, counts, cov


def test_separated_leading_eigenpairs_converge_and_come_back_descending():
    from icepack2_tools.laplace import arpack_action_bound, arpack_ncv
    top = [1000.0, 500.0, 250.0, 125.0]
    solve, counts, cov = _diagonal_problem(top + list(1.0 + 1e-3 * np.arange(300)))
    vals, vecs, converged = solve(4, tol=1e-10, maxiter=20)
    assert converged
    assert vals == pytest.approx(top, rel=1e-8)
    assert counts["H"] <= arpack_action_bound(4, arpack_ncv(4), 20)
    gram = vecs.T @ (vecs / cov[:, None])            # v_i . Gamma^-1 v_j
    assert gram == pytest.approx(np.eye(4), abs=1e-6)


def test_the_restart_limit_keeps_the_converged_pairs():
    """Two separated eigenvalues above a dense cluster: one restart cycle
    converges the two and cannot resolve the cluster to machine precision, so
    the solve returns those two instead of raising, within its action bound."""
    from icepack2_tools.laplace import arpack_action_bound
    solve, counts, cov = _diagonal_problem([100.0, 50.0] + list(1.0 + 1e-3 * np.arange(400)))
    vals, vecs, converged = solve(6, tol=0.0, maxiter=1, ncv=20)
    assert not converged
    assert 2 <= vals.size < 6
    assert vals[:2] == pytest.approx([100.0, 50.0], rel=1e-8)
    assert vecs.shape == (402, vals.size)
    assert counts["H"] <= arpack_action_bound(6, 20, 1)
    gram = vecs.T @ (vecs / cov[:, None])
    assert gram == pytest.approx(np.eye(vals.size), abs=1e-6)


def test_no_converged_pair_is_an_error():
    solve, counts, cov = _diagonal_problem(list(1.0 + 1e-4 * np.arange(400)))
    with pytest.raises(RuntimeError, match="no eigenpair"):
        solve(3, tol=0.0, maxiter=1, ncv=4)


def _dense_problem(n=9, m=11, seed=11):
    """A small linear-Gaussian problem with dense actions: full-rank data so
    the preconditioned operator's spectrum lies strictly above 1."""
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(n, n))
    prior_cov = A @ A.T + n * np.eye(n)
    prior_prec = np.linalg.inv(prior_cov)
    G = rng.normal(size=(m, n))
    W = np.diag(rng.uniform(0.5, 4.0, size=m))
    H = G.T @ W @ G
    g = rng.normal(size=n)
    exact = g @ np.linalg.solve(H + prior_prec, g)
    lam, V = scipy.linalg.eigh(H, prior_prec)               # V^T Gamma^-1 V = I
    order = np.argsort(-lam)
    counts = {"H": 0}

    def H_action(x):
        counts["H"] += 1
        return H @ x

    actions = (H_action, lambda x: prior_prec @ x, lambda r: prior_cov @ r)
    return dict(n=n, prior_cov=prior_cov, H=H, g=g, exact=exact, lam=lam[order],
                V=V[:, order], actions=actions, counts=counts)


def test_conjugate_gradients_bracket_the_variance_and_converge():
    from icepack2_tools.laplace import posterior_variance_cg
    P = _dense_problem()
    out = posterior_variance_cg(*P["actions"], P["g"], iterations=P["n"])
    lower, upper = out["lower"], out["upper"]
    assert out["iterations"] == P["n"] and P["counts"]["H"] == P["n"]
    assert np.all(np.diff(lower) >= -1e-9 * P["exact"])
    assert np.all(lower <= P["exact"] * (1 + 1e-9))
    assert np.all(upper >= P["exact"] * (1 - 1e-9))
    assert upper[0] == pytest.approx(P["g"] @ P["prior_cov"] @ P["g"])     # the prior variance
    assert lower[-1] == pytest.approx(P["exact"], rel=1e-8)
    assert upper[-1] == pytest.approx(P["exact"], rel=1e-6)
    assert P["g"] @ out["x"] == pytest.approx(P["exact"], rel=1e-8)


def test_the_eigen_seed_starts_at_the_truncated_variance_and_tightens_it():
    from icepack2_tools.laplace import posterior_variance, posterior_variance_cg
    P = _dense_problem()
    k = 3
    modes, lam = P["V"][:, :k].T, P["lam"][:k]
    truncated = posterior_variance(P["g"] @ P["prior_cov"] @ P["g"], modes @ P["g"], lam)
    out = posterior_variance_cg(*P["actions"], P["g"], iterations=P["n"] - k, modes=modes,
                                eigenvalues=lam)
    assert out["upper"][0] == pytest.approx(truncated, rel=1e-10)
    assert out["lower"][0] == pytest.approx(np.sum((modes @ P["g"]) ** 2 / (1 + lam)),
                                            rel=1e-10)
    assert out["projections"] == pytest.approx(modes @ P["g"])
    assert P["counts"]["H"] == P["n"] - k + 1                    # one action seeds
    assert out["lower"][-1] == pytest.approx(P["exact"], rel=1e-8)
    assert out["upper"][-1] == pytest.approx(P["exact"], rel=1e-6)
    assert np.all(out["lower"] <= P["exact"] * (1 + 1e-9))
    assert np.all(out["upper"] >= P["exact"] * (1 - 1e-9))


def test_a_resumed_solve_continues_the_same_bracket():
    from icepack2_tools.laplace import posterior_variance_cg
    P = _dense_problem()
    whole = posterior_variance_cg(*P["actions"], P["g"], iterations=5)
    first = posterior_variance_cg(*P["actions"], P["g"], iterations=2)
    rest = posterior_variance_cg(*P["actions"], P["g"], iterations=3, state=first)
    assert rest["iterations"] == 5 and rest["actions"] == 3
    assert rest["lower"] == pytest.approx(whole["lower"], rel=1e-10)
    assert rest["upper"] == pytest.approx(whole["upper"], rel=1e-8)


def test_stop_and_callback_follow_the_steps():
    from icepack2_tools.laplace import posterior_variance_cg
    P = _dense_problem()
    seen = []
    out = posterior_variance_cg(*P["actions"], P["g"], iterations=8, stop=lambda k: k >= 3,
                                callback=seen.append)
    assert out["iterations"] == 3 and [r["iteration"] for r in seen] == [0, 1, 2, 3]
    assert [r["lower"] for r in seen] == pytest.approx(list(out["lower"]))
    assert [r["upper"] for r in seen] == pytest.approx(list(out["upper"]))
    assert list(seen[-1]["lower_by_step"]) == pytest.approx(list(out["lower"]))
    whole = posterior_variance_cg(*P["actions"], P["g"], iterations=5)
    rest = posterior_variance_cg(*P["actions"], P["g"], iterations=3, state=seen[2])
    assert rest["iterations"] == 5
    assert rest["upper"] == pytest.approx(whole["upper"], rel=1e-8)


def test_a_prior_only_direction_keeps_the_bounds_valid():
    """Data that leave a direction unobserved put an eigenvalue of exactly 1
    in the preconditioned operator, the Gauss-Radau node itself; the bracket
    must survive the resulting near-singular quadrature matrices."""
    from icepack2_tools.laplace import posterior_variance_cg
    P = _dense_problem(n=8, m=5)
    out = posterior_variance_cg(*P["actions"], P["g"], iterations=P["n"])
    assert np.all(np.isfinite(out["upper"]))
    assert np.all(out["lower"] <= P["exact"] * (1 + 1e-9))
    assert np.all(out["upper"] >= P["exact"] * (1 - 1e-9))
    assert out["lower"][-1] == pytest.approx(P["exact"], rel=1e-8)


def test_quadrature_bounds_match_the_iterate_norm_and_refuse_bad_input():
    from icepack2_tools.laplace import cg_jacobi_matrix, posterior_variance_cg, quadrature_bounds
    P = _dense_problem()
    out = posterior_variance_cg(*P["actions"], P["g"], iterations=4)
    lower, upper = quadrature_bounds(out["alphas"], out["betas"], out["rho"][0])
    assert lower == pytest.approx(out["lower"]) and upper == pytest.approx(out["upper"])
    T = cg_jacobi_matrix(out["alphas"], out["betas"])
    ritz = np.linalg.eigvalsh(T)
    assert np.all(ritz > 1.0)                       # within the operator's spectrum
    with pytest.raises(ValueError, match="beta"):
        quadrature_bounds([1.0, 2.0], [0.5], 1.0)
    with pytest.raises(ValueError, match="positive"):
        quadrature_bounds([1.0], [0.5], 1.0, lower_eigenvalue=0.0)


fd = pytest.importorskip("firedrake")


def _slab(direction, weight_value=None):
    """A 4 x 1 strip of DG0 ice 100 m thick, grounded for x < 2, under a
    uniform velocity of 1 m/yr along `direction`."""
    from icepack2_tools.laplace import grounded_discharge_forms
    mesh = fd.RectangleMesh(4, 1, 4.0, 1.0)
    x, _ = fd.SpatialCoordinate(mesh)
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    h = fd.Function(Q0).assign(100.0)
    grounded = fd.Function(Q0).interpolate(fd.conditional(x < 2.0, 1.0, 0.0))
    u = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1)).assign(fd.Constant(direction))
    weight = None if weight_value is None else fd.Function(Q0).assign(weight_value)
    gl, edge = grounded_discharge_forms(u, h, grounded, weight)
    return fd.assemble(gl), fd.assemble(edge)


def test_the_grounding_line_flux_is_positive_from_grounded_to_floating():
    gl, edge = _slab((1.0, 0.0))
    assert gl == pytest.approx(100.0)        # h u across the unit facet at x = 2
    assert edge == pytest.approx(0.0)         # grounded ice leaves nowhere else


def test_flow_onto_grounded_ice_counts_negative_and_leaves_at_the_edge():
    gl, edge = _slab((-1.0, 0.0))
    assert gl == pytest.approx(-100.0)
    assert edge == pytest.approx(100.0)       # out of the grounded cells at x = 0


def test_a_weight_on_the_grounded_side_scales_its_discharge():
    gl, edge = _slab((1.0, 0.0), weight_value=2.0)
    assert gl == pytest.approx(200.0)


def _bed_slab(downstream_bed, downstream_thickness):
    """A 4 x 1 strip, marine grounded ice 100 m thick on a bed at -50 m for
    x < 2 and the given bed and thickness beyond, under 1 m/yr along x; the
    grounded mask is built from the geometry."""
    from icepack2_tools.laplace import grounded_discharge_forms, grounded_indicator
    rho_ratio = 917.0 / 1028.0
    mesh = fd.RectangleMesh(4, 1, 4.0, 1.0)
    x, _ = fd.SpatialCoordinate(mesh)
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    b = fd.Function(Q0).interpolate(fd.conditional(x < 2.0, -50.0, downstream_bed))
    h = fd.Function(Q0).interpolate(fd.conditional(x < 2.0, 100.0, downstream_thickness))
    s = fd.Function(Q0).interpolate(fd.max_value(b + h, (1 - rho_ratio) * h))
    grounded = grounded_indicator(Q0, s, b, rho_ratio)
    u = fd.Function(fd.VectorFunctionSpace(mesh, "CG", 1)).assign(fd.Constant((1.0, 0.0)))
    gl, edge = grounded_discharge_forms(u, h, grounded)
    return fd.assemble(gl), fd.assemble(edge)


def test_flow_into_a_floating_marine_cell_is_grounding_line_flux():
    gl, edge = _bed_slab(downstream_bed=-500.0, downstream_thickness=100.0)
    assert gl == pytest.approx(100.0)
    assert edge == pytest.approx(0.0)


def test_flow_onto_ice_free_land_books_nowhere():
    gl, edge = _bed_slab(downstream_bed=10.0, downstream_thickness=0.0)
    assert gl == pytest.approx(0.0)
    assert edge == pytest.approx(0.0)


def test_the_discharge_space_check_reads_the_element_not_a_context_key():
    """A run's context names its geometry space differently across checkouts
    (the production line's has no ``geometry_space`` key), so the check on the
    thickness space reads the element itself."""
    from icepack2_tools.laplace import is_dg0
    mesh = fd.UnitSquareMesh(2, 2)
    assert is_dg0(fd.FunctionSpace(mesh, "DG", 0))
    assert not is_dg0(fd.FunctionSpace(mesh, "DG", 1))
    assert not is_dg0(fd.FunctionSpace(mesh, "CG", 1))


def test_gauss_newton_actions_release_their_tangent_linears():
    """Every Gauss-Newton action makes a tangent-linear copy of each forward
    variable and keeps it in a dictionary on that variable until the map's
    finaliser removes it. tlm_adjoint 0.9.1dev (April 2026) registers a
    variable for that removal only when it creates the dictionary, so from
    the second action on the copies stay attached: one copy of the forward
    state per action, which ended a 100-mode eigensolve out of memory. Our
    tlm_adjoint (github.com/hoffmaao/tlm_adjoint, branch ismip7) registers
    the variable on every write; this test fails on an environment that
    still installs upstream."""
    import gc
    tla = pytest.importorskip("tlm_adjoint.firedrake")

    mesh = fd.UnitSquareMesh(4, 4)
    V = fd.FunctionSpace(mesh, "CG", 1)
    x, y = fd.SpatialCoordinate(mesh)
    m = fd.Function(V, name="m").interpolate(0.1 * x)
    f = fd.Function(V, name="f").interpolate(fd.sin(fd.pi * x) * fd.sin(fd.pi * y))
    u = fd.Function(V, name="u")
    v = fd.TestFunction(V)
    F = fd.inner(fd.grad(u), fd.grad(v)) * fd.dx + fd.exp(m) * u * v * fd.dx - f * v * fd.dx
    bc = fd.DirichletBC(V, 0.0, "on_boundary")
    tla.reset_manager()
    tla.start_manager()
    tla.EquationSolver(F == 0, u, bc).solve()
    tla.stop_manager()
    H = tla.CachedGaussNewton(u, lambda tau: fd.assemble(fd.inner(tau, v) * fd.dx))
    dm = fd.Function(V, name="dm").assign(1.0)

    def attached():
        return sum(len(getattr(z, "_tlm_adjoint__tangent_linears", {})) for z in (m, f, u))

    results = []
    for _ in range(3):
        ddJ, = H.action([m], [dm])
        results.append(ddJ.dat.data_ro.copy())
        del ddJ
        gc.collect()
        assert attached() == 0
    assert np.allclose(results[1], results[0]) and np.allclose(results[2], results[0])
    assert np.linalg.norm(results[0]) > 0.0


def _save_modes(path, mesh, values):
    Q = fd.FunctionSpace(mesh, "CG", 1)
    x, y = fd.SpatialCoordinate(mesh)
    with fd.CheckpointFile(str(path), "w") as chk:
        chk.save_mesh(mesh)
        for i, (ca, cp) in enumerate(values):
            chk.save_function(fd.Function(Q, name=f"mode_{i:04d}_alpha").interpolate(ca * x + y))
            chk.save_function(fd.Function(Q, name=f"mode_{i:04d}_phi").interpolate(cp * y - x))
        chk.set_attr("/", "n_modes", len(values))
    return Q


def test_saved_modes_come_back_on_the_live_mesh(tmp_path):
    """The modes load onto the mesh that saved them by a dof copy, onto the
    same mesh under another numbering by interpolation, and refuse another
    mesh altogether."""
    from icepack2_tools.laplace import load_modes
    path = tmp_path / "eigendec.h5"
    mesh = fd.UnitSquareMesh(4, 4)
    Q = _save_modes(path, mesh, [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)])
    modes, copied = load_modes(str(path), Q, count=2)
    assert copied and len(modes) == 2
    x, y = fd.SpatialCoordinate(mesh)
    want = fd.Function(Q).interpolate(3.0 * x + y)
    assert np.allclose(modes[1][0].dat.data_ro, want.dat.data_ro)
    assert modes[1][1].name() == "mode_0001_phi"

    other = fd.UnitSquareMesh(4, 4, reorder=False)
    Q_other = fd.FunctionSpace(other, "CG", 1)
    modes, copied = load_modes(str(path), Q_other)
    assert not copied and len(modes) == 3
    xo, yo = fd.SpatialCoordinate(other)
    want = fd.Function(Q_other).interpolate(6.0 * yo - xo)
    assert np.allclose(modes[2][1].dat.data_ro, want.dat.data_ro, atol=1e-10)

    with pytest.raises(ValueError, match="another mesh"):
        load_modes(str(path), fd.FunctionSpace(fd.UnitSquareMesh(5, 5), "CG", 1))
