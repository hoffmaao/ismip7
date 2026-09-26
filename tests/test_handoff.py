r"""One objective across restarts (icepack2_tools.handoff): the settings
comparison, the accepted-iterate lookup, and the objective handoff gap.
Pure Python, no mesh."""
import numpy as np
import pytest

from icepack2_tools.handoff import (
    OBJECTIVE_KEYS,
    accepted_evaluation,
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
        "prior_sigma_alpha": 0.1828, "prior_rho_theta": 7500.0,
        "friction_control": "sqrt", "friction": "budd", "n_flow": 3.0,
        "geometry_space": "dg0", "friction_anchor_length": 20000.0,
        "lake_ice_base": 1, "fluidity_prior_origin": "thermomechanical",
        "grad_precond": "mass_consistent",
        "subelement_friction": 0, "exact_front": 0,
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
