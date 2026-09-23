r"""check_budd_map's forward re-solve compares the forward's velocity with the
MAP's on two mesh objects that Firedrake numbers differently. A raw dat copy
compared permuted fields and reported rel L2 0.57 for a forward that matched
the MAP to 1e-7 (Sep 20 2026). The comparison aligns nodes by coordinate."""
import importlib.util
import os

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _check_budd_map():
    spec = importlib.util.spec_from_file_location(
        "check_budd_map", os.path.join(REPO, "antarctica", "scripts", "check_budd_map.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_a_permuted_numbering_is_undone():
    m = _check_budd_map()
    rng = np.random.default_rng(1)
    x = rng.uniform(-3e6, 3e6, size=(500, 2))
    field = rng.normal(size=(500, 2))
    perm = rng.permutation(500)
    ia, ib = m.node_permutation(x, x[perm])
    recovered = np.empty_like(field)
    recovered[ia] = field[perm][ib]
    assert np.array_equal(recovered, field)
    # the defect: the same data taken in the other mesh's order
    assert not np.array_equal(field[perm], field)


def test_coordinates_that_round_the_same_still_match():
    m = _check_budd_map()
    x = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    ia, ib = m.node_permutation(x, x[::-1] + 1e-9)
    assert np.allclose(x[ia], (x[::-1] + 1e-9)[ib])


def test_a_different_vertex_set_is_refused():
    m = _check_budd_map()
    x = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    with pytest.raises(SystemExit):
        m.node_permutation(x, x + 0.5)
    with pytest.raises(SystemExit):
        m.node_permutation(x, x[:2])
