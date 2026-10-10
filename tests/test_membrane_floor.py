r"""The membrane floor (``h_visc_floor``): the inversion's knob, and the floor
a forward takes from its MAP (``runconfig.forward_hvisc_floor``). At 10 m the
first water row beside the ice carried the ocean drag onto the front (issue
#153), so the inversion runs 2.5 m and records it; every MAP from before the
record was inverted at 10 m."""
import pytest

from icepack2_tools.runconfig import forward_hvisc_floor, hvisc_floor


def test_the_inversion_runs_two_and_a_half_metres_unless_told(monkeypatch):
    monkeypatch.delenv("ISMIP7_RC_HVISC_FLOOR", raising=False)
    assert hvisc_floor() == 2.5
    monkeypatch.setenv("ISMIP7_RC_HVISC_FLOOR", "10")
    assert hvisc_floor() == 10.0


def test_a_forward_runs_the_floor_its_map_records(monkeypatch):
    monkeypatch.delenv("ISMIP7_RC_HVISC_FLOOR", raising=False)
    assert forward_hvisc_floor(1.0) == 1.0
    assert forward_hvisc_floor("10.0") == 10.0
    monkeypatch.setenv("ISMIP7_RC_HVISC_FLOOR", "1")
    assert forward_hvisc_floor(1.0) == 1.0
    with pytest.raises(RuntimeError, match="follows its MAP"):
        forward_hvisc_floor(10.0, source="m.h5")


def test_a_map_without_the_record_runs_ten_metres_or_the_knob(monkeypatch):
    monkeypatch.delenv("ISMIP7_RC_HVISC_FLOOR", raising=False)
    assert forward_hvisc_floor(None) == 10.0
    monkeypatch.setenv("ISMIP7_RC_HVISC_FLOOR", "1")
    assert forward_hvisc_floor(None) == 1.0
