r"""Two run knobs that every driver reads through icepack2_tools.runconfig.

``ISMIP7_FRICTION`` is a closed set: theta means a different thing under each
law, and a mistyped value used to run the legacy action branch without a word.
``ISMIP7_MESH=checkpoint`` names the mesh embedded in the MAP or restart file,
which is the only way a job that sources site_env.sh (submit.sh projection)
can run MAP-native, since that file always exports a derived path.
"""
import pytest

from icepack2_tools.runconfig import (
    FRICTION_LAWS, MESH_FROM_CHECKPOINT, friction, mesh_override,
)


@pytest.mark.parametrize("law", FRICTION_LAWS)
def test_every_law_in_the_set_is_accepted(monkeypatch, law):
    monkeypatch.setenv("ISMIP7_FRICTION", law)
    assert friction() == law


def test_spelling_is_normalised_before_the_check(monkeypatch):
    monkeypatch.setenv("ISMIP7_FRICTION", " Budd ")
    assert friction() == "budd"


@pytest.mark.parametrize("bad", ["rc", "coulomb", "weertman", ""])
def test_an_unknown_law_is_an_error_at_startup(monkeypatch, bad):
    monkeypatch.setenv("ISMIP7_FRICTION", bad)
    with pytest.raises(ValueError, match="ISMIP7_FRICTION must be one of"):
        friction()


def test_the_default_law_is_in_the_set(monkeypatch):
    monkeypatch.delenv("ISMIP7_FRICTION", raising=False)
    assert friction() in FRICTION_LAWS


@pytest.mark.parametrize("value", [None, "", "  ", MESH_FROM_CHECKPOINT])
def test_unset_empty_and_the_sentinel_all_mean_the_checkpoint_mesh(monkeypatch, value):
    if value is None:
        monkeypatch.delenv("ISMIP7_MESH", raising=False)
    else:
        monkeypatch.setenv("ISMIP7_MESH", value)
    assert mesh_override() is None


def test_a_path_is_the_compute_mesh(monkeypatch):
    monkeypatch.setenv("ISMIP7_MESH", "/meshes/antarctica_10000_1000_buffered20000.msh")
    assert mesh_override() == "/meshes/antarctica_10000_1000_buffered20000.msh"
