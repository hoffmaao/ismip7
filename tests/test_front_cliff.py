r"""The push on an ice front inside the mesh against the terminus condition
on a front that is the mesh boundary (issue #153).

ISMIP7's cell-wise residual on a uniform slab of DG0 ice in x < L: one mesh
ends at x = L with the right side a calving id, the other runs on with
ice-free cells. The velocity block evaluated at the zero state against
v = e_x is the net front force (every interior jump vanishes on a uniform
slab). The facet driving stress alone matches the terminus condition for
floating ice and on flat land, falls short at a grounded marine cliff by
g D (rho_I H - rho_W D) / 2, and beside rock of another height is set by the
step in bed height; exact_front gives every such edge the terminus push.
"""
import numpy as np
import pytest

fd = pytest.importorskip("firedrake")

from finat.ufl import FiniteElement  # noqa: E402
from icepack2.constants import gravity as g, ice_density as rho_I, water_density as rho_W  # noqa: E402
from icepack2_tools.dual_friction import build_rc_residual  # noqa: E402

L, B, W, NX = 20e3, 10e3, 4e3, 10


def _front_force(H_ice, b_val, buffered, exact_front, h_water=0.0, cg1_geometry=False,
                 b_water=None):
    mesh = (fd.RectangleMesh(int(NX * (L + B) / L), 2, L + B, W) if buffered
            else fd.RectangleMesh(NX, 2, L, W))
    x, _y = fd.SpatialCoordinate(mesh)
    dg0 = FiniteElement("DG", "triangle", 0)
    Z = (fd.VectorFunctionSpace(mesh, "CG", 1)
         * fd.TensorFunctionSpace(mesh, dg0, symmetry=True)
         * fd.VectorFunctionSpace(mesh, dg0))
    z = fd.Function(Z)
    Qg = fd.FunctionSpace(mesh, "CG", 1) if cg1_geometry else fd.FunctionSpace(mesh, "DG", 0)
    H = fd.Function(Qg).interpolate(fd.conditional(fd.lt(x, L), H_ice, h_water))
    b = fd.Function(Qg).interpolate(
        fd.conditional(fd.lt(x, L), b_val, b_val if b_water is None else b_water))
    s = fd.Function(Qg).interpolate(fd.max_value(b + H, (1 - rho_I / rho_W) * H))
    Q = fd.FunctionSpace(mesh, "CG", 1)
    one = fd.Function(Q).assign(1.0)
    F = build_rc_residual(
        z, fd.Function(Q), fd.Function(Q), H=H, s=s, b=b, C_w0=one, A4_base=one,
        n_flow=fd.Constant(3.0), n_flow_val=3.0, m_slide=3.0, tau_c=0.1, alpha=1e-4,
        H_ref=100.0, k_lim=0.0, calving_ids=(2,), exact_front=exact_front, front_hmin=1.0)
    r = fd.assemble(F).subfunctions[0].dat.data_ro
    return r[:, 0].sum() / W


def _exact(H_ice, b_val):
    s = max(b_val + H_ice, (1 - rho_I / rho_W) * H_ice)
    d = max(0.0, H_ice - s)
    return 0.5 * g * (rho_I * H_ice ** 2 - rho_W * d ** 2)


CASES = [("floating", 400.0, -1000.0), ("grounded marine cliff", 1600.0, -300.0),
         ("deep grounded marine cliff", 800.0, -600.0), ("land", 1000.0, 100.0)]


@pytest.mark.parametrize("name, H_ice, b_val", CASES)
def test_the_terminus_condition_is_the_exact_push(name, H_ice, b_val):
    assert _front_force(H_ice, b_val, buffered=False, exact_front=False) == \
        pytest.approx(_exact(H_ice, b_val), rel=1e-12)


@pytest.mark.parametrize("name, H_ice, b_val", CASES)
def test_the_facet_push_falls_short_only_at_a_grounded_marine_cliff(name, H_ice, b_val):
    internal = _front_force(H_ice, b_val, buffered=True, exact_front=False)
    exact = _exact(H_ice, b_val)
    grounded_marine = b_val < 0 and b_val + H_ice > (1 - rho_I / rho_W) * H_ice
    gap = 0.5 * g * (-b_val) * (rho_I * H_ice + rho_W * b_val) if grounded_marine else 0.0
    assert exact - internal == pytest.approx(gap, rel=1e-9, abs=1e-9 * exact)
    if grounded_marine:
        assert internal < 0.9 * exact


@pytest.mark.parametrize("name, H_ice, b_val", CASES)
def test_exact_front_gives_an_internal_cliff_the_terminus_push(name, H_ice, b_val):
    assert _front_force(H_ice, b_val, buffered=True, exact_front=True) == \
        pytest.approx(_exact(H_ice, b_val), rel=1e-12)


def test_a_film_on_the_water_side_keeps_the_correction_exact():
    r"""Water cells holding 0.4 m, below ISMIP7_FRONT_HMIN: still open water to
    the correction, and the push it completes is the facet term's own, with
    the film in avg(H)."""
    H_ice, b_val = 1600.0, -300.0
    pushed = _front_force(H_ice, b_val, buffered=True, exact_front=True, h_water=0.4)
    # the film's own facet with the boundary beyond pushes too (0.4 m, zero
    # in practice); everything else is the exact cliff push
    assert pushed == pytest.approx(_exact(H_ice, b_val), rel=1e-6)


@pytest.mark.parametrize("b_rock", [0.0, 300.0, 700.0])
def test_a_land_margin_beside_rock_of_another_height_gets_the_cliff_push(b_rock):
    r"""Ice 500 m thick on a 100 m bed (surface at 600 m) beside bare rock at
    0 m (the facet push 1.2 times the cliff push outward), at 300 m (0.6 times
    outward) and at 700 m, above the ice surface (0.2 times back into the
    ice). exact_front replaces the bed-step push with the free-cliff push in
    each case; restricting it to ocean facets kept the bed-step push and
    raised the misfit of IU's final Budd MAP at its own controls by 3.5 %
    (issue #153). Against rock above the ice surface the cliff push is a known
    error (issue #166)."""
    H_ice, b_ice = 500.0, 100.0
    facet = _front_force(H_ice, b_ice, buffered=True, exact_front=False, b_water=b_rock)
    pushed = _front_force(H_ice, b_ice, buffered=True, exact_front=True, b_water=b_rock)
    assert pushed == pytest.approx(_exact(H_ice, b_ice), rel=1e-6)
    assert abs(facet - _exact(H_ice, b_ice)) > 0.1 * _exact(H_ice, b_ice)


def test_exact_front_refuses_cg1_geometry():
    with pytest.raises(ValueError, match="DG0"):
        _front_force(400.0, -1000.0, buffered=True, exact_front=True, cg1_geometry=True)
