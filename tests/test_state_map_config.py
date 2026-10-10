r"""A forward state checkpoint carries the MAP's control and grounding-scheme
attributes (``simulation.save_model_state``), so a chained link restarted from
it rebuilds the residual the MAP was inverted under."""
import importlib.util
import os
import sys

import firedrake as fd
import pytest
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


@pytest.mark.parametrize("exact_front", [1, 2])
def test_state_carries_map_configuration(tmp_path, exact_front):
    sim = _simulation()
    metadata = {"friction_control": "exp", "friction_c_ref": 0.25,
                "subelement_friction": 1, "exact_front": exact_front,
                "fluidity_control": "floating", "drag_gate": "vertex",
                "h_visc_floor": 1.0}
    path = str(tmp_path / "state.h5")
    sim.save_model_state(_ctx(metadata), path, 1.0)
    got = _attrs(path, sim.MAP_CONFIG_KEYS)
    assert got["friction_control"] == "exp"
    assert float(got["friction_c_ref"]) == 0.25
    assert int(got["subelement_friction"]) == 1
    assert int(got["exact_front"]) == exact_front
    assert got["fluidity_control"] == "floating"
    assert got["drag_gate"] == "vertex"
    assert float(got["h_visc_floor"]) == 1.0


def test_state_omits_what_the_map_did_not_record(tmp_path):
    sim = _simulation()
    path = str(tmp_path / "state.h5")
    sim.save_model_state(_ctx({}), path, 1.0)
    assert _attrs(path, sim.MAP_CONFIG_KEYS) == {}


def test_state_carries_the_initial_state_and_the_relaxation_record(tmp_path):
    r"""Every forward checkpoint names the initial state its chain began from
    and, for a relaxed MAP, how the MAP's geometry was made
    (icepack2_tools.relaxation), so a restart and a report know it."""
    sim = _simulation()
    ctx = _ctx({})
    ctx["init_state"] = "relaxed-controls"
    ctx["relax_record"] = {"relax_t_start": 2014.0, "relax_forcing": "ocx protocol",
                           "relax_state": "relax_x_2000_final.h5"}
    path = str(tmp_path / "state.h5")
    sim.save_model_state(ctx, path, 2003.0)
    got = _attrs(path, ("init_state", "relax_t_start", "relax_forcing", "relax_state",
                        "relaxation_end_state"))
    assert got["init_state"] == "relaxed-controls"
    assert float(got["relax_t_start"]) == 2014.0
    assert got["relax_forcing"] == "ocx protocol"
    assert got["relax_state"] == "relax_x_2000_final.h5"
    # only the relaxation driver's own final state marks itself an end state
    assert "relaxation_end_state" not in got
    final = str(tmp_path / "final.h5")
    sim.save_model_state(_ctx({}), final, 2015.0,
                         extra_attrs={"stalled": 0, "relaxation_end_state": 1})
    assert int(_attrs(final, ("relaxation_end_state",))["relaxation_end_state"]) == 1
