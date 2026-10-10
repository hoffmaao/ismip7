r"""The push on an ice front inside the mesh against the terminus condition
on a front that is the mesh boundary (issue #153).

ISMIP7's cell-wise residual on a uniform slab of DG0 ice in x < L: one mesh
ends at x = L with the right side a calving id, the other runs on with
ice-free cells. The velocity block evaluated at the zero state against
v = e_x is the net front force (every interior jump vanishes on a uniform
slab). The facet driving stress alone matches the terminus condition for
floating ice and on flat land, falls short at a grounded marine cliff by
g D (rho_I H - rho_W D) / 2, and beside rock of another height is set by the
step in bed height; exact_front gives every such edge the terminus push
(version 1), or the push of the face above the ice-free neighbour's bed
(version 2, issue #166), which leaves none against rock above the ice surface.
"""
import pytest

fd = pytest.importorskip("firedrake")

from finat.ufl import FiniteElement  # noqa: E402
from icepack2.constants import gravity as g, ice_density as rho_I, water_density as rho_W  # noqa: E402
from icepack2_tools.dual_friction import build_rc_residual, front_cliff_correction  # noqa: E402
from icepack2_tools.runconfig import exact_front_version, forward_exact_front  # noqa: E402

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


@pytest.mark.parametrize("version", [1, 2])
def test_a_film_on_the_water_side_keeps_the_correction_exact(version):
    r"""Water cells holding 0.4 m, below ISMIP7_FRONT_HMIN: still open water to
    the correction, and the push it completes is the facet term's own, with
    the film in avg(H)."""
    H_ice, b_val = 1600.0, -300.0
    pushed = _front_force(H_ice, b_val, buffered=True, exact_front=version, h_water=0.4)
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
    (issue #153). Against rock above the ice surface version 1's cliff push
    drives the ice into the rock; version 2 below leaves none (issue #166)."""
    H_ice, b_ice = 500.0, 100.0
    facet = _front_force(H_ice, b_ice, buffered=True, exact_front=False, b_water=b_rock)
    pushed = _front_force(H_ice, b_ice, buffered=True, exact_front=True, b_water=b_rock)
    assert pushed == pytest.approx(_exact(H_ice, b_ice), rel=1e-6)
    assert abs(facet - _exact(H_ice, b_ice)) > 0.1 * _exact(H_ice, b_ice)


@pytest.mark.parametrize("version", [True, 2])
def test_exact_front_refuses_cg1_geometry(version):
    with pytest.raises(ValueError, match="DG0"):
        _front_force(400.0, -1000.0, buffered=True, exact_front=version, cg1_geometry=True)


def _exposed(H_ice, b_ice, b_rock):
    r"""Version 2's push: the face above the neighbour's bed ``b_rock``."""
    s = max(b_ice + H_ice, (1 - rho_I / rho_W) * H_ice)
    d = max(0.0, H_ice - s)
    h_e = min(H_ice, max(s - b_rock, 0.0))
    d_e = min(d, max(-b_rock, 0.0))
    return 0.5 * g * (rho_I * h_e ** 2 - rho_W * d_e ** 2)


@pytest.mark.parametrize("name, H_ice, b_val", CASES)
def test_version_2_keeps_the_free_cliff_push_beside_a_bed_at_the_ice_base(name, H_ice, b_val):
    r"""Floating ice over deep water and a cliff whose neighbour's bed is its own:
    versions 1 and 2 give the same push."""
    assert _front_force(H_ice, b_val, buffered=True, exact_front=2) == \
        pytest.approx(_exact(H_ice, b_val), rel=1e-12)


@pytest.mark.parametrize("b_rock", [-50.0, 0.0, 100.0, 101.0, 300.0, 599.0, 600.0, 700.0, 2000.0])
def test_version_2_pushes_with_the_face_above_the_neighbour_bed(b_rock):
    r"""500 m of ice on a 100 m bed (surface at 600 m) beside bare ground from
    below the ice base to far above the surface: the full push down to the
    base, the exposed face between base and surface, and none against rock at
    or above the surface, continuously across both ends."""
    H_ice, b_ice = 500.0, 100.0
    pushed = _front_force(H_ice, b_ice, buffered=True, exact_front=2, b_water=b_rock)
    assert pushed == pytest.approx(_exposed(H_ice, b_ice, b_rock), rel=1e-6,
                                   abs=1e-9 * _exact(H_ice, b_ice))
    if b_rock >= 600.0:
        assert abs(pushed) < 1e-9 * _exact(H_ice, b_ice)


def test_version_2_worked_example():
    r"""The numbers the docstring and the README quote: 0.36 of the push beside
    rock at 300 m, all of it at 0 m, none at 700 m."""
    H_ice, b_ice = 500.0, 100.0
    full = _exact(H_ice, b_ice)
    assert _exposed(H_ice, b_ice, 0.0) / full == pytest.approx(1.0)
    assert _exposed(H_ice, b_ice, 300.0) / full == pytest.approx(0.36)
    assert _exposed(H_ice, b_ice, 700.0) == 0.0


@pytest.mark.parametrize("H_ice, b_ice, b_shoal", [(400.0, -1000.0, -200.0), (1600.0, -300.0, -100.0)])
def test_version_2_leaves_water_only_above_a_shoal(H_ice, b_ice, b_shoal):
    r"""A seabed shoaler than the ice base beside floating ice and beside a
    grounded marine cliff: the face below the shoal rests on it, and the water
    stands from sea level down to the shoal."""
    pushed = _front_force(H_ice, b_ice, buffered=True, exact_front=2, b_water=b_shoal)
    assert pushed == pytest.approx(_exposed(H_ice, b_ice, b_shoal), rel=1e-6)
    assert 0.0 < pushed < _exact(H_ice, b_ice)


def test_version_2_needs_the_bed():
    mesh = fd.UnitSquareMesh(2, 2)
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    v = fd.TestFunction(fd.VectorFunctionSpace(mesh, "CG", 1))
    with pytest.raises(ValueError, match="needs the bed"):
        front_cliff_correction(v, fd.Function(Q0), fd.Function(Q0), mesh, version=2)


@pytest.mark.parametrize("value, version", [
    (False, 0), (True, 1), (0, 0), (1, 1), (2, 2), (2.0, 2), ("2", 2), (" 1 ", 1), (b"2", 2)])
def test_exact_front_versions(value, version):
    assert exact_front_version(value) == version


@pytest.mark.parametrize("value", [3, -1, 1.5, "on", "", None])
def test_exact_front_refuses_an_unknown_version(value):
    with pytest.raises(ValueError, match="exact_front"):
        exact_front_version(value)


def test_a_forward_follows_the_maps_cliff_push(monkeypatch):
    monkeypatch.delenv("ISMIP7_EXACT_FRONT", raising=False)
    assert forward_exact_front(2) == 2
    assert forward_exact_front(None) == 0              # a MAP from before the record
    monkeypatch.setenv("ISMIP7_EXACT_FRONT", "2")
    assert forward_exact_front(2) == 2                 # the knob may repeat it
    with pytest.raises(RuntimeError, match="follows its MAP"):
        forward_exact_front(1)                         # and may not change it
    with pytest.raises(RuntimeError, match="follows its MAP"):
        forward_exact_front(None)
