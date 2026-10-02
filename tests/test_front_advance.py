r"""A retreat-only front (ISMIP7_FRONT_ADVANCE=none): the law retreats the
front and nothing advances past the t=0 extent."""
import numpy as np
import pytest

from icepack2_tools import runconfig
from icepack2_tools.front import front_removal_mask


def test_the_knob_parses_and_refuses_the_rest(monkeypatch):
    monkeypatch.delenv("ISMIP7_FRONT_ADVANCE", raising=False)
    assert runconfig.front_advance() == "free"
    monkeypatch.setenv("ISMIP7_FRONT_ADVANCE", "None")
    assert runconfig.front_advance() == "none"
    monkeypatch.setenv("ISMIP7_FRONT_ADVANCE", "retreat")
    with pytest.raises(ValueError):
        runconfig.front_advance()


def test_retreat_only_removes_beyond_the_initial_extent_as_well():
    law = np.array([False, True, False, False])     # the front passed cell 1
    init = np.array([False, False, True, True])     # cells 2, 3 were ice-free at t=0
    assert np.array_equal(front_removal_mask(law, init, False), law)
    assert np.array_equal(front_removal_mask(law, init, True),
                          np.array([False, True, True, True]))
    assert np.array_equal(front_removal_mask(law, None, True), law)
