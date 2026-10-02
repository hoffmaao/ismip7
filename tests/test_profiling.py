r"""The span accumulator behind the inversion's timing record (issue 156)."""
import pytest

pytest.importorskip("mpi4py")
from mpi4py import MPI  # noqa: E402

from icepack2_tools.profiling import Spans  # noqa: E402


def test_a_span_entered_twice_reports_its_total_and_reduce_clears():
    spans = Spans(MPI.COMM_WORLD)
    spans.add("prior_solve", 1.5)
    spans.add("gather", 0.25)
    spans.add("prior_solve", 2.0)
    with spans("residual_norm"):
        pass
    out = spans.reduce()
    assert list(out) == ["gather", "prior_solve", "residual_norm"]
    assert out["prior_solve"] == pytest.approx(3.5)
    assert out["gather"] == pytest.approx(0.25)
    assert 0.0 <= out["residual_norm"] < 1.0
    assert spans.reduce() == {}


def test_a_span_is_recorded_when_its_block_raises():
    spans = Spans(MPI.COMM_WORLD)
    with pytest.raises(ValueError):
        with spans("forward"):
            raise ValueError("diverged")
    assert "forward" in spans.reduce()


def test_clear_drops_the_spans_without_reducing():
    spans = Spans(MPI.COMM_WORLD)
    spans.add("set_controls", 1.0)
    spans.clear()
    assert spans.reduce() == {}
