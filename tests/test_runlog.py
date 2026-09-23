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
    ({"id": "x", "task": "inversion", "title": "t", "status": "done"}, "owner"),
    ({"id": "x", "task": "nonsense", "title": "t", "status": "done",
      "owner": "o"}, "task"),
    ({"id": "x", "task": "inversion", "title": "t", "status": "finished",
      "owner": "o"}, "status"),
    ({"id": "x", "task": "inversion", "title": "t", "status": "done",
      "owner": "o", "resolution": "2 km"}, "unknown field"),
    ({"id": "other", "task": "inversion", "title": "t", "status": "done",
      "owner": "o"}, "filename"),
])
def test_a_malformed_record_is_refused(tmp_path, bad, message):
    m = _build_runlog()
    (tmp_path / "x.json").write_text(json.dumps(bad))
    with pytest.raises(ValueError, match=message):
        m.load(tmp_path)


def test_the_csv_carries_every_field(tmp_path):
    m = _build_runlog()
    out = tmp_path / "log.csv"
    n = m.write_csv(m.load(), out)
    import csv
    rows = list(csv.reader(open(out)))
    assert rows[0] == [h for _, h in m.FIELDS]
    assert len(rows) == n + 1
    # a list field joins into one cell rather than breaking the row
    assert all(len(r) == len(m.FIELDS) for r in rows)
