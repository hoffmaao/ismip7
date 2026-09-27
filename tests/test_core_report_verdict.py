r"""The core report states a check it could not run as unjudged.

check_ismip6_track.py and compare_ismip6.py exit 2 when they cannot judge a
run, and compare_ismip6.py does so on every machine without the ISMIP6
ensemble archive (neither Quartz nor the Mac holds one). The report read
every non-zero status as a failure, so an unjudged run went on record as off
track or outside the envelope.
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

import core_report  # noqa: E402


def test_a_pass_and_a_failure_keep_their_words():
    assert core_report.verdict(0, "ON TRACK", "OFF TRACK") == "ON TRACK"
    assert core_report.verdict(1, "ON TRACK", "OFF TRACK") == "OFF TRACK"
    assert core_report.verdict(0, "inside envelope", "outside envelope") == "inside envelope"
    assert core_report.verdict(1, "inside envelope", "outside envelope") == "outside envelope"


def test_a_check_that_could_not_run_is_neither():
    words = core_report.verdict(2, "inside envelope", "outside envelope")
    assert "outside" not in words and "inside" not in words
    assert words.startswith("not judged") and "exit 2" in words
