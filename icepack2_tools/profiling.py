r"""Wall-clock spans inside one inversion evaluation, reduced across ranks.

The inversion's timing record splits an evaluation into the forward, the
adjoint and the rest (issue #156). ``Spans`` names the rest: each span is timed
with ``perf_counter`` on every rank and reported as the SLOWEST rank's time,
so the number printed from rank 0 describes the communicator's critical path
(AGENTS.md section 3). Accompanying durations and an unspanned remainder use
the same convention. The remainder is formed from rank-local measurements
before its maximum is reduced. The reduction is collective; call
:meth:`Spans.reduce` on every rank at the same point of the evaluation.
"""

from contextlib import contextmanager
from time import perf_counter

import numpy as np

from icepack2_tools.mpi_stats import global_max


class Spans:
    r"""Accumulates named wall-clock spans on this rank until reduced.

    A name entered twice before a reduce accumulates, so a span that runs once
    per control or once per line-search trial reports its total.
    """

    def __init__(self, comm):
        self.comm = comm
        self._seconds = {}

    @contextmanager
    def __call__(self, name):
        t0 = perf_counter()
        try:
            yield
        finally:
            self.add(name, perf_counter() - t0)

    def add(self, name, seconds):
        self._seconds[name] = self._seconds.get(name, 0.0) + float(seconds)

    def clear(self):
        r"""Drop this rank's spans; local, no communication."""
        self._seconds.clear()

    def reduce(self, *, durations=None, unspanned_total=None):
        r"""``{name: slowest rank's seconds}`` in sorted name order, then
        clear. Collective. Every rank must hold the same names, which holds
        whenever the code the spans wrap runs the same on every rank; a rank
        that disagrees raises on every rank rather than pairing one rank's
        span with another's.

        ``durations`` supplies rank-local durations reduced by the same
        maximum. ``unspanned_total`` supplies a rank-local duration containing
        every accumulated span; its remainder is computed locally, then
        returned as ``unspanned`` among the reduced durations. Calls using
        either option return ``(spans, durations)``.
        """
        durations_given = durations is not None or unspanned_total is not None
        local_durations = {
            name: float(seconds) for name, seconds in (durations or {}).items()
        }
        if unspanned_total is not None:
            if "unspanned" in local_durations:
                raise ValueError("unspanned is reserved for unspanned_total")
            local_durations["unspanned"] = (
                float(unspanned_total) - sum(self._seconds.values())
            )

        span_names = sorted(self._seconds)
        duration_names = sorted(local_durations)
        signature = (tuple(span_names), tuple(duration_names))
        if len(set(self.comm.allgather(signature))) != 1:
            raise RuntimeError(
                "Spans.reduce: the ranks supplied different timing names; "
                "the timed code did not run the same on every rank")
        out = {
            name: global_max(np.array([self._seconds[name]]), comm=self.comm)
            for name in span_names
        }
        reduced_durations = {
            name: global_max(np.array([local_durations[name]]), comm=self.comm)
            for name in duration_names
        }
        self._seconds.clear()
        return (out, reduced_durations) if durations_given else out
