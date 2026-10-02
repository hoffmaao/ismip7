r"""A forward state checkpoint carries the MAP's control and grounding-scheme
attributes (``simulation.save_model_state``), so a chained link restarted from
it rebuilds the residual the MAP was inverted under."""
import importlib.util
import os
import sys

import firedrake as fd
from firedrake import (Function, FunctionSpace, TensorFunctionSpace,
                       UnitSquareMesh, VectorFunctionSpace, FiniteElement)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "antarctica", "scripts")


def _simulation():
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    spec = importlib.util.spec_from_file_location(
        "simulation", os.path.join(SCRIPTS, "simulation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _ctx(metadata):
    mesh = UnitSquareMesh(2, 2)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Z = V * TensorFunctionSpace(mesh, dg0, symmetry=True) * VectorFunctionSpace(mesh, dg0)
    Q = FunctionSpace(mesh, "CG", 1)
    Q0 = FunctionSpace(mesh, "DG", 0)
    return {
        "mesh": mesh, "z": Function(Z), "h": Function(Q0),
        "theta": Function(Q), "phi": Function(Q), "u_obs": Function(V),
        "b": Function(Q0), "s": Function(Q0), "phi_eff": Function(Q0),
        "geom_dg": True, "checkpoint_metadata": metadata,
    }


def _attrs(path, keys):
    with fd.CheckpointFile(path, "r") as chk:
        return {k: chk.get_attr("/", k) for k in keys if chk.has_attr("/", k)}


def test_state_carries_map_configuration(tmp_path):
    sim = _simulation()
    metadata = {"friction_control": "exp", "friction_c_ref": 0.25,
                "subelement_friction": 1, "exact_front": 1,
                "fluidity_control": "floating"}
    path = str(tmp_path / "state.h5")
    sim.save_model_state(_ctx(metadata), path, 1.0)
    got = _attrs(path, sim.MAP_CONFIG_KEYS)
    assert got["friction_control"] == "exp"
    assert float(got["friction_c_ref"]) == 0.25
    assert int(got["subelement_friction"]) == 1
    assert int(got["exact_front"]) == 1
    assert got["fluidity_control"] == "floating"


def test_state_omits_what_the_map_did_not_record(tmp_path):
    sim = _simulation()
    path = str(tmp_path / "state.h5")
    sim.save_model_state(_ctx({}), path, 1.0)
    assert _attrs(path, sim.MAP_CONFIG_KEYS) == {}
