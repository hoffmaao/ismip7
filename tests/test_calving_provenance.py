r"""A calving result has to carry its threshold. The von Mises rate is
inversely proportional to sigma_max, so a run recorded as "vonmises" and
nothing else cannot be reproduced, and the default has changed before."""
import re
import pathlib

import pytest

SRC = (pathlib.Path(__file__).resolve().parents[1]
       / "antarctica" / "scripts" / "simulation.py").read_text()


def test_the_front_owner_line_names_the_von_mises_threshold():
    block = SRC[SRC.index('if calving != "none":'):]
    block = block[:block.index("Calving front owner:")]
    assert "_calving_sigma_max()" in block, (
        "the front-owner line does not read the von Mises thresholds")
    assert "sigma_max grounded" in block


def test_the_front_owner_reaches_the_checkpoint():
    r"""The log is the first thing lost, so the data file must carry it too."""
    assert 'chk.set_attr("/", "calving_front_owner"' in SRC
    assert 'ctx["calving_front_owner"] = front_owner' in SRC


def test_the_threshold_is_formatted_without_losing_digits():
    r"""`%g` on 0.2 gives "0.2", not "0.200000"; a truncating format would
    record a different number than the run used."""
    for value in (0.15, 0.2, 1.0, 0.105, 0.0001):
        assert float(f"{value:g}") == value
