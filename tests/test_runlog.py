r"""The tracked run log. Every simulation leaves a record in
``antarctica/runlog/``, and ``reports/SIMULATIONS.md`` is rendered from them,
so a record that does not validate or a render that has gone stale is a test
failure rather than something a reader discovers later."""
import importlib.util
import json
import os

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNLOG = os.path.join(REPO, "antarctica", "runlog")


def _build_runlog():
    spec = importlib.util.spec_from_file_location(
        "build_runlog",
        os.path.join(REPO, "antarctica", "scripts", "build_runlog.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_every_record_validates():
    m = _build_runlog()
    records = m.load()
    assert records, "the run log is empty"
    ids = [r["id"] for r in records]
    assert len(ids) == len(set(ids))


def test_the_rendered_log_is_current():
    m = _build_runlog()
    assert m.main(["--check"]) == 0


@pytest.mark.parametrize("bad,message", [
    ({"id": "x", "task": "inversion", "title": "t", "status": "done"}, "institution"),
    ({"id": "x", "task": "nonsense", "title": "t", "status": "done",
      "institution": "o"}, "task"),
    ({"id": "x", "task": "inversion", "title": "t", "status": "finished",
      "institution": "o"}, "status"),
    ({"id": "x", "task": "inversion", "title": "t", "status": "done",
      "institution": "o", "resolution": "2 km"}, "unknown field"),
    ({"id": "other", "task": "inversion", "title": "t", "status": "done",
      "institution": "o"}, "filename"),
])
def test_a_malformed_record_is_refused(tmp_path, bad, message):
    m = _build_runlog()
    (tmp_path / "x.json").write_text(json.dumps(bad))
    with pytest.raises(ValueError, match=message):
        m.load(tmp_path)


@pytest.mark.parametrize("payload", ["[]", "3", '"a string"'])
def test_a_record_that_is_not_an_object_is_refused(tmp_path, capsys, payload):
    m = _build_runlog()
    (tmp_path / "x.json").write_text(payload)
    assert m.main(["--runlog", str(tmp_path), "--check",
                   "--output", str(tmp_path / "out.md")]) == 2
    assert "x.json: not a JSON object" in capsys.readouterr().err


def test_a_rewrite_with_no_record_change_leaves_the_log_unchanged(tmp_path):
    m = _build_runlog()
    out = tmp_path / "SIMULATIONS.md"
    assert m.main(["--write", "--output", str(out)]) == 0
    first = out.read_text()
    assert m.main(["--write", "--output", str(out)]) == 0
    assert out.read_text() == first
