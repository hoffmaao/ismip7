r"""``mpi_stats.global_rss_mib``: the per-evaluation memory figure the
inversion writes into its timing record (issue 159)."""
import pytest

MPI = pytest.importorskip("mpi4py.MPI")

from icepack2_tools.mpi_stats import global_rss_mib  # noqa: E402


def test_rss_over_one_rank():
    rss = global_rss_mib(MPI.COMM_WORLD)
    assert set(rss) == {"mean", "max", "peak_max"}
    assert rss["mean"] is not None and rss["mean"] > 1.0
    assert rss["max"] == pytest.approx(rss["mean"])
    # the peak is the high-water mark the current RSS has reached
    assert rss["peak_max"] >= 0.99 * rss["max"]


def test_rss_sees_memory_python_does_not_track():
    import numpy as np

    before = global_rss_mib(MPI.COMM_WORLD)["mean"]
    block = np.ones(64 * 2**20 // 8)  # 64 MiB, written so it is resident
    after = global_rss_mib(MPI.COMM_WORLD)["mean"]
    assert after - before > 48.0
    del block
