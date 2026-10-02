r"""One objective across restarts (icepack2_tools.handoff): the settings
comparison, the accepted-iterate lookup, and the objective handoff gap.
Pure Python, no mesh."""
import numpy as np
import pytest

from icepack2_tools.handoff import (
    OBJECTIVE_KEYS,
    accepted_evaluation,
    frozen_in_control,
    handoff_gap,
    objective_mismatches,
)


def _settings(**over):
    base = {
        "misfit_norm": "sigma", "misfit_scale": 918512.0,
        "log_vel_weight": 85380.44865839917, "log_vel_eps": 1.0,
        "dhdt_weight": 1.0, "dhdt_net_sigma": 10.0,
        "prior_form": "bilaplacian", "gamma_theta": 1e5, "gamma_phi": 1e5,
        "prior_sigma_theta": 0.3, "prior_sigma_phi": 0.3, "prior_rho": 7500.0,
        "prior_sigma_alpha": 0.1828, "prior_rho_theta": 7500.0, "friction_c_ref": 0.0334,
        "friction_control": "sqrt", "friction": "budd", "n_flow": 3.0,
        "geometry_space": "dg0", "friction_anchor_length": 20000.0,
        "lake_ice_base": 1, "fluidity_prior_origin": "thermomechanical",
        "grad_precond": "mass_consistent",
        "subelement_friction": 1, "subelement_scheme": "sep1", "exact_front": 0,
        "fluidity_control": "all",
    }
    base.update(over)
    return base


def test_identical_settings_have_no_mismatch():
    assert objective_mismatches(_settings(), _settings()) == []


def test_every_objective_setting_is_compared():
    for key in OBJECTIVE_KEYS:
        current = _settings()
        current[key] = "changed" if isinstance(current[key], str) else current[key] * 2 + 1
        out = objective_mismatches(_settings(), current)
        assert len(out) == 1 and out[0].startswith(f"{key}:"), key


def test_numbers_are_compared_as_numbers_and_missing_keys_are_not_mismatches():
    recorded = _settings(prior_rho="7500.0")          # an attribute read back as text
    assert objective_mismatches(recorded, _settings()) == []
    old = _settings()
    del old["misfit_scale"], old["prior_sigma_alpha"]  # an older checkpoint
    assert objective_mismatches(old, _settings()) == []
    # the auto weight re-derived on a restart is exactly the drift to catch
    out = objective_mismatches(_settings(), _settings(log_vel_weight=56646.35))
    assert out == ["log_vel_weight: 85380.44865839917 -> 56646.35"]


def test_a_sub_element_record_without_a_scheme_was_sep2():
    legacy = _settings(subelement_friction=1)
    del legacy["subelement_scheme"]
    out = objective_mismatches(legacy, _settings(subelement_friction=1, subelement_scheme="sep1"))
    assert out == ["subelement_scheme: 'sep2' -> 'sep1'"]
    assert objective_mismatches(legacy, _settings(subelement_friction=1,
                                                  subelement_scheme="sep2")) == []
    # without sub-element friction the scheme selects nothing
    off = dict(subelement_friction=0)
    plain = _settings(**off)
    del plain["subelement_scheme"]
    assert objective_mismatches(plain, _settings(**off)) == []
    assert objective_mismatches(_settings(subelement_scheme="sep2", **off),
                                _settings(subelement_scheme="sep1", **off)) == []
    no_record = _settings(subelement_scheme="sep2")
    del no_record["subelement_friction"]
    assert objective_mismatches(no_record, _settings(subelement_scheme="sep1", **off)) == []


def test_the_accepted_evaluation_is_the_one_matching_the_objective():
    ring = [(1.0e9, "trial a"), (7.5e8, "accepted"), (7.9e8, "trial b")]
    assert accepted_evaluation(ring, 7.5e8)[1] == "accepted"
    assert accepted_evaluation(ring, 7.5e8 * (1 + 1e-12))[1] == "accepted"
    assert accepted_evaluation(ring, 7.0e8) is None
    assert accepted_evaluation([], 1.0) is None


def test_the_handoff_gap_is_relative():
    assert handoff_gap(1.0e9, 1.0e9) == 0.0
    assert handoff_gap(1.0e9, 1.001e9) == pytest.approx(1e-3)
    assert np.isfinite(handoff_gap(0.0, 1.0))


def test_auto_weight_frozen_only_within_one_control():
    sqrt_record = _settings()
    assert frozen_in_control(sqrt_record, "prior_sigma_alpha", "sqrt") == 0.1828
    assert frozen_in_control(sqrt_record, "prior_sigma_alpha", "exp") is None
    exp_record = _settings(friction_control="exp", prior_sigma_alpha=1.0)
    assert frozen_in_control(exp_record, "prior_sigma_alpha", "sqrt") is None
    assert frozen_in_control(exp_record, "friction_c_ref", "exp") == 0.0334
    legacy = {"prior_sigma_alpha": 0.2}
    assert frozen_in_control(legacy, "prior_sigma_alpha", "sqrt") is None
    assert frozen_in_control(_settings(friction_c_ref=0.0), "friction_c_ref", "sqrt") is None
