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
