r"""A posterior-variance resume continues only under the likelihood its saved
vectors were computed under (run_posterior_variance.check_resume); state saved
with no likelihood record counts as the velocity term alone."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "antarctica", "scripts"))
# the script pulls icepack2 -> irksome, which must be imported before any UFL form is assembled
import run_posterior_variance as rpv  # noqa: E402

import json  # noqa: E402

import numpy as np  # noqa: E402

VELOCITY = {"likelihood": "velocity", "terms": ["velocity"]}
MAP = {"likelihood": "map", "terms": ["velocity", "dhdt"], "dhdt_weight": 1.0}


def _saved_state(path, likelihood=None):
    extra = {} if likelihood is None else {"likelihood": json.dumps(likelihood)}
    np.savez(path, n_modes=0, eigendec="", eigenvalues=np.zeros(0), lower_eigenvalue=1.0, **extra)
    return np.load(path)


def _check(z, likelihood):
    return rpv.check_resume(z, 0, "", np.zeros(0), likelihood, "state.npz")


def test_a_resume_under_another_likelihood_is_refused(tmp_path):
    with _saved_state(tmp_path / "v.npz", VELOCITY) as z:
        assert _check(z, VELOCITY) is None
        assert "likelihood" in _check(z, MAP)
    with _saved_state(tmp_path / "m.npz", MAP) as z:
        assert _check(z, MAP) is None
        assert "likelihood" in _check(z, VELOCITY)
        assert "likelihood" in _check(z, dict(MAP, dhdt_weight=2.0))


def test_state_saved_with_no_likelihood_resumes_under_the_velocity_term_alone(tmp_path):
    with _saved_state(tmp_path / "old.npz") as z:
        assert _check(z, VELOCITY) is None
        assert "likelihood" in _check(z, MAP)


def test_a_resume_with_other_modes_is_still_refused(tmp_path):
    with _saved_state(tmp_path / "v.npz", VELOCITY) as z:
        assert "modes" in rpv.check_resume(z, 2, "eigendec.h5", np.ones(2), VELOCITY, "state.npz")
