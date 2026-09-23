r"""A diagnostic solve that converges onto a speed no ice can have must be
treated as a failed solve, so the step reaches the rescue ladder and the
subcycles instead of being carried forward.

Measured motivation: a 1 km control in the Lambert/Amery grounding trough
reported SNES reason=2 with a function norm below 1 and a peak speed of
2.5e6 m/yr, and the step loop accepted it, because it asked only whether the
solver had converged.
"""
import pathlib

import pytest

from icepack2_tools.runconfig import max_speed_bound

SRC = (pathlib.Path(__file__).resolve().parents[1]
       / "antarctica" / "scripts" / "simulation.py").read_text()


def test_the_default_bound_is_above_any_real_ice_and_below_a_runaway(monkeypatch):
    monkeypatch.delenv("ISMIP7_MAX_SPEED", raising=False)
    bound = max_speed_bound()
    assert bound == pytest.approx(20000.0)
    # the fastest ice measured anywhere is about 17 km/yr, and a healthy 1 km
    # control peaks near 5.5, so the bound must sit above both
    assert bound > 17000.0
    # and below the speeds the trough reached once it ignited
    assert bound < 20435.0 * 1.001


def test_the_bound_can_be_disabled_but_not_inverted(monkeypatch):
    monkeypatch.setenv("ISMIP7_MAX_SPEED", "0")
    assert max_speed_bound() == 0.0
    monkeypatch.setenv("ISMIP7_MAX_SPEED", "-1")
    with pytest.raises(ValueError, match="positive"):
        max_speed_bound()
    monkeypatch.setenv("ISMIP7_MAX_SPEED", "5e4")
    assert max_speed_bound() == pytest.approx(5e4)


def test_a_disabled_bound_short_circuits_before_interpolating():
    r"""The check runs after every accepted solve, so with the bound off it
    must not pay for an interpolation and a global reduction."""
    body = SRC[SRC.index("def _speed_bound_violation"):]
    body = body[:body.index("def _record_transport")]
    assert body.index("if bound <= 0.0:") < body.index("speed_diag.interpolate")


def test_the_direct_path_rejects_an_unphysical_solve():
    block = SRC[SRC.index('solve_diagnostic(f"step-{k}-direct")'):]
    block = block[:block.index("except fd.ConvergenceError as exc:")]
    assert "_speed_bound_violation" in block
    assert "raise fd.ConvergenceError" in block, (
        "a violation on the direct path must enter the rescue ladder")


def test_every_rescue_rung_is_checked_too():
    r"""A rung that converges onto the same runaway must not be accepted just
    because it was reached by trust region."""
    block = SRC[SRC.index("                try:\n                    action()"):]
    block = block[:block.index("except fd.ConvergenceError:")]
    assert "_speed_bound_violation" in block
    assert block.index("_speed_bound_violation") < block.index("return True")


def test_the_violation_is_reported_with_its_location():
    body = SRC[SRC.index("def _speed_bound_violation"):]
    body = body[:body.index("def _record_transport")]
    assert "ISMIP7_MAX_SPEED" in body, "the message must name the knob"
    assert "u_max_xy" in body, "a runaway is local; say where it is"
