r"""The horizontal-force-balance calving criterion.

``icepack2_tools.calving_laws`` is the criterion of Buck (2023), Coffey et al.
(2024), Coffey and Lai (2025) and Slater and Wagner (2025): a front holds while
the resistive stress it carries is below what a crevasse field can bear. The
numbers below are the ones those papers pin, so a change to the algebra shows
up here rather than as a front that quietly stops calving.
"""
import pytest
from firedrake import (Constant, Function, FunctionSpace, TensorFunctionSpace,
                       UnitSquareMesh, VectorFunctionSpace, FiniteElement,
                       as_matrix, as_vector)
from icepack2.constants import gravity as G, ice_density as RHO_I, water_density as RHO_W

from icepack2_tools.calving_laws import (buttressing_number,
                                         critical_stress,
                                         density_in_model_units,
                                         height_above_flotation,
                                         hfb_calving_rate, resistive_stress)

H_SHELF = 300.0
#: a bed deep enough that 300 m of ice floats freely
B_DEEP = -2000.0


@pytest.fixture(scope="module")
def spaces():
    mesh = UnitSquareMesh(4, 4)
    dg0 = FiniteElement("DG", "triangle", 0)
    return (mesh, FunctionSpace(mesh, "DG", 0),
            TensorFunctionSpace(mesh, dg0, symmetry=True),
            VectorFunctionSpace(mesh, "DG", 0))


def _cells(space, expr):
    return Function(space).interpolate(expr).dat.data_ro


def test_the_density_conversion_is_icepack2s_own():
    assert density_in_model_units(1024.0) == pytest.approx(float(RHO_W), rel=1e-12)
    assert density_in_model_units(917.0) == pytest.approx(float(RHO_I), rel=1e-12)


def test_a_freely_floating_shelf_has_zero_height_above_flotation(spaces):
    mesh, Q0, _, _ = spaces
    # h chosen so the ice floats exactly: h = (rho_w / rho_i) * (-b)
    b = Constant(-100.0)
    h = Constant(float(RHO_W / RHO_I) * 100.0)
    assert _cells(Q0, height_above_flotation(h, b)).max() == pytest.approx(0.0, abs=1e-9)


def test_the_resistive_stress_is_the_front_normal_component(spaces):
    r"""``n . M n``: with the front normal along x it is M_xx, which is the
    papers' ``R_xx``, and not the larger principal value."""
    mesh, Q0, Sigma, W0 = spaces
    M = Function(Sigma).interpolate(as_matrix(((0.2, 0.05), (0.05, 0.4))))
    n = Function(W0).interpolate(as_vector((1.0, 0.0)))
    assert _cells(Q0, resistive_stress(M, n)).max() == pytest.approx(0.2)
    n_y = Function(W0).interpolate(as_vector((0.0, 1.0)))
    assert _cells(Q0, resistive_stress(M, n_y)).max() == pytest.approx(0.4)


def test_an_unbuttressed_shelf_sits_exactly_at_the_threshold(spaces):
    r"""Buck's result, and Coffey and Lai's ``B = 0``: at zero tensile
    strength the critical stress equals the depth-averaged extensional stress
    of a freely spreading shelf, ``rho_i g H (1 - rho_i/rho_w) / 2``. That
    equality is the whole criterion, so it is the first thing to check."""
    mesh, Q0, _, _ = spaces
    h, b = Constant(H_SHELF), Constant(B_DEEP)
    r_crit = _cells(Q0, critical_stress(h, b, sigma_max=0.0))
    expected = (float(RHO_I) * float(G) * H_SHELF
                * (1.0 - float(RHO_I / RHO_W)) / 2.0)
    assert r_crit.max() == pytest.approx(expected, rel=1e-10)


def test_tensile_strength_raises_the_threshold(spaces):
    mesh, Q0, _, _ = spaces
    h, b = Constant(H_SHELF), Constant(B_DEEP)
    weak = _cells(Q0, critical_stress(h, b, 0.0)).max()
    strong = _cells(Q0, critical_stress(h, b, 0.15)).max()
    assert strong > weak


def test_fresh_water_in_the_crevasse_lowers_the_threshold(spaces):
    r"""A meltwater-filled basal crevasse is easier to open than a seawater
    one, so the front is closer to failure."""
    mesh, Q0, _, _ = spaces
    h, b = Constant(H_SHELF), Constant(B_DEEP)
    sea = _cells(Q0, critical_stress(h, b, 0.0)).max()
    fresh = _cells(Q0, critical_stress(
        h, b, 0.0, rho_c=density_in_model_units(1000.0))).max()
    assert fresh < sea


def test_buttressing_is_zero_at_failure_and_one_under_no_stress(spaces):
    mesh, Q0, Sigma, W0 = spaces
    h, b = Constant(H_SHELF), Constant(B_DEEP)
    n = Function(W0).interpolate(as_vector((1.0, 0.0)))
    free = float(_cells(Q0, critical_stress(h, b, 0.0)).max())
    at_failure = Function(Sigma).interpolate(as_matrix(((free, 0.0), (0.0, 0.0))))
    assert _cells(Q0, buttressing_number(at_failure, n, h, b)).max() \
        == pytest.approx(0.0, abs=1e-9)
    none = Function(Sigma).interpolate(as_matrix(((0.0, 0.0), (0.0, 0.0))))
    assert _cells(Q0, buttressing_number(none, n, h, b)).max() \
        == pytest.approx(1.0, rel=1e-12)


def test_a_fully_buttressed_front_does_not_calve(spaces):
    r"""The threshold behaviour, which is what separates this law from a
    proportional one: below failure it removes nothing at all."""
    mesh, Q0, Sigma, W0 = spaces
    h, b = Constant(H_SHELF), Constant(B_DEEP)
    n = Function(W0).interpolate(as_vector((1.0, 0.0)))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((1000.0, 0.0)))
    M = Function(Sigma).interpolate(as_matrix(((0.0, 0.0), (0.0, 0.0))))
    assert _cells(Q0, hfb_calving_rate(u, M, h, b, n)).max() \
        == pytest.approx(0.0, abs=1e-12)


def test_at_failure_the_front_retreats_at_the_ice_speed(spaces):
    r"""``c = |u| (R_xx/R_crit)^1``, so a front exactly at the threshold
    retreats at the speed the ice arrives, which is a stationary front."""
    mesh, Q0, Sigma, W0 = spaces
    h, b = Constant(H_SHELF), Constant(B_DEEP)
    n = Function(W0).interpolate(as_vector((1.0, 0.0)))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((1000.0, 0.0)))
    free = float(_cells(Q0, critical_stress(h, b, 0.0)).max())
    M = Function(Sigma).interpolate(as_matrix(((free, 0.0), (0.0, 0.0))))
    assert _cells(Q0, hfb_calving_rate(u, M, h, b, n)).max() \
        == pytest.approx(1000.0, rel=1e-6)


def test_the_ratio_cap_bounds_what_one_step_removes(spaces):
    r"""``R_crit`` scales with the thickness, so on ice thinning to nothing the
    ratio is unbounded; without the cap one step could carry the front across
    many cells."""
    mesh, Q0, Sigma, W0 = spaces
    h, b = Constant(H_SHELF), Constant(B_DEEP)
    n = Function(W0).interpolate(as_vector((1.0, 0.0)))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((100.0, 0.0)))
    M = Function(Sigma).interpolate(as_matrix(((50.0, 0.0), (0.0, 0.0))))
    assert _cells(Q0, hfb_calving_rate(u, M, h, b, n, ratio_max=5.0)).max() \
        == pytest.approx(500.0, rel=1e-6)


def test_compression_does_not_calve(spaces):
    mesh, Q0, Sigma, W0 = spaces
    h, b = Constant(H_SHELF), Constant(B_DEEP)
    n = Function(W0).interpolate(as_vector((1.0, 0.0)))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((500.0, 0.0)))
    M = Function(Sigma).interpolate(as_matrix(((-0.3, 0.0), (0.0, 0.0))))
    assert _cells(Q0, hfb_calving_rate(u, M, h, b, n)).max() \
        == pytest.approx(0.0, abs=1e-12)


def test_the_thickness_rule_holds_the_front_at_its_critical_thickness(spaces):
    r"""``Hc`` is where the front settles, not a cutoff: at ``H = Hc`` the rate
    is exactly the speed the ice arrives with, so the front is stationary. That
    is what makes the rule self-limiting where a stress threshold is not."""
    from icepack2_tools.calving_laws import thickness_calving_rate
    mesh, Q0, _, _ = spaces
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((500.0, 0.0)))
    hc, b = 150.0, Constant(B_DEEP)
    assert _cells(Q0, thickness_calving_rate(u, Constant(hc), b, hc)).max() \
        == pytest.approx(500.0, rel=1e-6)
    # thinner than Hc: faster than the ice arrives, so the front retreats
    assert _cells(Q0, thickness_calving_rate(u, Constant(50.0), b, hc)).max() > 500.0
    # thicker: slower, so the front advances
    thick = _cells(Q0, thickness_calving_rate(u, Constant(300.0), b, hc)).max()
    assert 0.0 <= thick < 500.0


def test_the_thickness_rule_never_calves_backwards(spaces):
    r"""Well above ``Hc`` the factor would go negative, which would read as the
    front advancing under its own calving law rather than under the ice flux."""
    from icepack2_tools.calving_laws import thickness_calving_rate
    mesh, Q0, _, _ = spaces
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((500.0, 0.0)))
    assert _cells(Q0, thickness_calving_rate(
        u, Constant(4000.0), Constant(B_DEEP), 150.0)).max() \
        == pytest.approx(0.0, abs=1e-12)


def test_the_thickness_rule_spares_a_margin_on_a_bed_above_sea_level(spaces):
    r"""The level set anchors on every ice edge, so a land margin or nunatak
    thinner than ``Hc`` would erode and be booked as calving. The gate is the
    bed below sea level: a grounded marine cliff still calves."""
    from icepack2_tools.calving_laws import thickness_calving_rate
    mesh, Q0, _, _ = spaces
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((500.0, 0.0)))
    h = Constant(50.0)
    land = _cells(Q0, thickness_calving_rate(u, h, Constant(200.0), 150.0))
    assert land.max() == pytest.approx(0.0, abs=1e-12)
    # 50 m of ice on a 10 m deep bed is grounded and still calves
    cliff = _cells(Q0, thickness_calving_rate(u, h, Constant(-10.0), 150.0))
    assert cliff.min() > 500.0


def test_the_critical_thickness_knob_rejects_a_nonpositive_value(monkeypatch):
    from icepack2_tools.runconfig import calving_thickness_hc
    monkeypatch.setenv("ISMIP7_CALVING_HC", "150")
    assert calving_thickness_hc() == pytest.approx(150.0)
    monkeypatch.setenv("ISMIP7_CALVING_HC", "0")
    with pytest.raises(ValueError, match="positive thickness"):
        calving_thickness_hc()


def test_the_floating_von_mises_default_is_the_published_antarctic_value(monkeypatch):
    r"""Wilner et al. (2023) calibrate this same law against ten Antarctic
    shelves and report 105-400 kPa, mean 225, median 230: the only one of the
    four laws they test whose parameter is consistent between shelves. The
    default is the centre of that cluster, so a silent change to it changes
    what the submission claims to be running."""
    from icepack2_tools.runconfig import calving_sigma_max
    monkeypatch.delenv("ISMIP7_CALVING_SIGMA_MAX_FLOATING", raising=False)
    monkeypatch.delenv("ISMIP7_CALVING_SIGMA_MAX_GROUNDED", raising=False)
    grounded, floating = calving_sigma_max()
    assert floating == pytest.approx(0.2)
    assert 0.105 <= floating <= 0.400, "outside the published Antarctic range"
    assert grounded == pytest.approx(1.0)


def test_the_critical_thickness_default_is_the_observed_front(monkeypatch):
    r"""The published Antarctic minimum-thickness thresholds are for a position
    law and are shelf-specific, so the transferable number for our rate form is
    the observed front thickness: BedMachine v4.1 has a 144.7 m mean and a
    152.4 m floating-front median."""
    from icepack2_tools.runconfig import calving_thickness_hc
    monkeypatch.delenv("ISMIP7_CALVING_HC", raising=False)
    assert calving_thickness_hc() == pytest.approx(150.0)


def test_the_resolved_calving_knobs_name_every_parameter_at_its_value(monkeypatch):
    r"""What the front owner line, the checkpoint and the core report record.
    A defaulted knob is absent from the environment, so it has to be stated
    at its value or two runs under different defaults read the same."""
    from icepack2_tools.runconfig import calving_knobs
    for k in ("SIGMA_MAX", "RHO_C", "HFB_EXPONENT", "HFB_RATIO_MAX", "HC",
              "SIGMA_MAX_GROUNDED", "SIGMA_MAX_FLOATING"):
        monkeypatch.delenv(f"ISMIP7_CALVING_{k}", raising=False)
    assert calving_knobs("hfb") == {
        "ISMIP7_CALVING_SIGMA_MAX": "0", "ISMIP7_CALVING_RHO_C": "1024",
        "ISMIP7_CALVING_HFB_EXPONENT": "1", "ISMIP7_CALVING_HFB_RATIO_MAX": "5"}
    assert calving_knobs("thickness") == {"ISMIP7_CALVING_HC": "150"}
    assert calving_knobs("vonmises") == {
        "ISMIP7_CALVING_SIGMA_MAX_GROUNDED": "1",
        "ISMIP7_CALVING_SIGMA_MAX_FLOATING": "0.2"}
    assert calving_knobs("none") == {}
    monkeypatch.setenv("ISMIP7_CALVING_RHO_C", "1000")
    assert calving_knobs("hfb")["ISMIP7_CALVING_RHO_C"] == "1000"
