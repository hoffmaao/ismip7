r"""A 2 km MAP can warm-start an inversion on the 1 km mesh: the helpers the
inversion uses to decide what the warm start supplies."""
import numpy as np
import pytest

firedrake = pytest.importorskip("firedrake")
from firedrake import Function, FunctionSpace, UnitSquareMesh  # noqa: E402

from icepack2_tools.transfer import interpolate_with_fill, meshes_match  # noqa: E402


def test_the_same_mesh_matches_and_a_finer_one_does_not():
    coarse, fine = UnitSquareMesh(4, 4), UnitSquareMesh(8, 8)
    assert meshes_match(coarse, coarse)
    assert meshes_match(coarse, UnitSquareMesh(4, 4))
    assert not meshes_match(coarse, fine)


def test_a_warm_start_from_a_smaller_mesh_fills_the_prior_with_cold_ice():
    r"""Target dofs the source mesh does not cover take the stated fill (the
    inversion passes the forward's cold-ice baseline for the fluidity prior),
    so A = A_prior exp(phi) stays positive there."""
    source_mesh = UnitSquareMesh(4, 4)
    target_mesh = firedrake.RectangleMesh(8, 8, 1.5, 1.5)     # a wider domain
    src = Function(FunctionSpace(source_mesh, "CG", 1)).assign(40.0)
    tgt = Function(FunctionSpace(target_mesh, "CG", 1))
    n_missing, n_total, n_clamped = interpolate_with_fill(tgt, src, 1.0)
    assert 0 < n_missing < n_total and n_clamped == 0
    vals = tgt.dat.data_ro
    assert set(np.round(vals, 9)) == {1.0, 40.0}
    assert vals.min() > 0.0


def test_a_map_s_own_mesh_loads_from_the_checkpoint_with_its_recorded_name(tmp_path):
    r"""Rice's 2 km MAPs travel without their .msh, so ISMIP7_MESH=checkpoint
    solves on the mesh the MAP carries, under the basename it recorded, which
    names the boundary-id sidecar."""
    from icepack2_tools.transfer import load_checkpoint_mesh
    mesh = UnitSquareMesh(4, 4, name="firedrake_default")
    path = str(tmp_path / "map.h5")
    with firedrake.CheckpointFile(path, "w") as chk:
        chk.save_mesh(mesh)
        chk.set_attr("/", "mesh_basename", "antarctica_5000_2000_buffered0.msh")
    loaded, basename = load_checkpoint_mesh(path)
    assert basename == "antarctica_5000_2000_buffered0.msh"
    assert meshes_match(loaded, mesh)


def test_a_checkpoint_that_names_no_mesh_is_refused(tmp_path):
    from icepack2_tools.transfer import load_checkpoint_mesh
    path = str(tmp_path / "map.h5")
    with firedrake.CheckpointFile(path, "w") as chk:
        chk.save_mesh(UnitSquareMesh(2, 2, name="firedrake_default"))
    with pytest.raises(ValueError, match="mesh_basename"):
        load_checkpoint_mesh(path)


def test_the_mesh_sentinel_reads_the_warm_start(monkeypatch):
    from icepack2_tools.runconfig import inversion_mesh_source
    monkeypatch.delenv("ISMIP7_MESH", raising=False)
    assert inversion_mesh_source("derived.msh") == ("derived.msh", False)
    monkeypatch.setenv("ISMIP7_MESH", "other.msh")
    assert inversion_mesh_source("derived.msh") == ("other.msh", False)
    monkeypatch.setenv("ISMIP7_MESH", "checkpoint")
    monkeypatch.setenv("ISMIP7_WARM_START", "rice_2km.h5")
    assert inversion_mesh_source("derived.msh") == ("rice_2km.h5", True)
    monkeypatch.delenv("ISMIP7_WARM_START")
    with pytest.raises(ValueError, match="ISMIP7_WARM_START"):
        inversion_mesh_source("derived.msh")


@pytest.mark.parametrize("value, expected",
                         [(None, True), ("", True), ("1", True), ("0", False)])
def test_eval_continuation_is_on_unless_zero(monkeypatch, value, expected):
    from icepack2_tools.runconfig import eval_continuation
    if value is None:
        monkeypatch.delenv("ISMIP7_EVAL_CONTINUATION", raising=False)
    else:
        monkeypatch.setenv("ISMIP7_EVAL_CONTINUATION", value)
    assert eval_continuation() is expected


def test_eval_continuation_refuses_a_word(monkeypatch):
    from icepack2_tools.runconfig import eval_continuation
    monkeypatch.setenv("ISMIP7_EVAL_CONTINUATION", "off")
    with pytest.raises(ValueError):
        eval_continuation()
