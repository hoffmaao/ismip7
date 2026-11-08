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
