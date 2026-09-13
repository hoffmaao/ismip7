r"""The adaptation counter lives in the mesh filename.

``adapt_mesh.py`` names each adapted mesh ``<root>_adapt<N>.msh`` and derives
the next name from the reference mesh it was handed.  It used to take ``N``
from a checkpoint ``adapt_count`` attribute, which only some writers stamp:
with the attribute absent it computed ``N = 0 + 1`` while ``root`` had already
lost the suffix, so ``X_adapt1`` regenerated its own name and the run either
overwrote the mesh it was reading or aborted on the guard.

Serial, no firedrake, no data files.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "antarctica", "scripts"))

from mesh_naming import adapt_lineage  # noqa: E402

# The two Ua-preset meshes and one parametric mesh actually in antarctica/mesh.
FRESH = [
    "antarctica_ua_180000_2000",
    "antarctica_ua_180000_2000_obs",
    "antarctica_64000_2500_buffered",
]


def next_mesh(basename):
    r"""What adapt_mesh.py writes when no --out-mesh is given."""
    root, k = adapt_lineage(basename)
    return f"{root}_adapt{k + 1}"


@pytest.mark.parametrize("name", FRESH)
def test_an_unadapted_mesh_starts_the_lineage(name):
    assert adapt_lineage(name) == (name, 0)
    assert next_mesh(name) == f"{name}_adapt1"


@pytest.mark.parametrize("name", FRESH)
def test_an_adapted_mesh_advances_instead_of_repeating(name):
    r"""The regression: a checkpoint carrying only mesh_basename (the MAP the
    inversion writes, and the forward's final.h5) must still advance."""
    adapted = f"{name}_adapt1"
    assert adapt_lineage(adapted) == (name, 1)
    assert next_mesh(adapted) == f"{name}_adapt2"
    assert next_mesh(adapted) != adapted


def test_the_lineage_never_revisits_a_name():
    r"""What the overwrite guard in adapt_mesh.py relies on: repeated
    adaptation is strictly increasing, so the output is never the reference."""
    name = FRESH[0]
    seen = [name]
    for _ in range(5):
        nxt = next_mesh(seen[-1])
        assert nxt not in seen
        seen.append(nxt)
    assert seen[-1] == f"{name}_adapt5"


def test_accepts_a_path_or_a_msh_filename():
    r"""adapt_mesh.py passes a bare basename; the inversion passes the full
    ISMIP7_MESH path it loaded."""
    name = FRESH[0]
    for spelling in (f"{name}_adapt3.msh",
                     f"/some/where/antarctica/mesh/{name}_adapt3.msh"):
        assert adapt_lineage(spelling) == (name, 3)


def test_only_a_trailing_counter_is_a_lineage_suffix():
    r"""`_adapt` inside a name is not a counter, and neither is a
    non-numeric suffix: both belong to the root."""
    assert adapt_lineage("antarctica_adapted_test_2000") == (
        "antarctica_adapted_test_2000", 0)
    assert adapt_lineage("antarctica_2000_adapt") == (
        "antarctica_2000_adapt", 0)
    assert adapt_lineage("antarctica_adapt2_2000") == (
        "antarctica_adapt2_2000", 0)
