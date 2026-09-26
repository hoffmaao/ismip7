r"""The exponent continuation a momentum solve is started from.

A single solve at the physical exponents can fail from a cold guess, so the
flow and sliding exponents are ramped from 1 to their targets in equal
steps, each step solved from the previous one. On a fine mesh with a rough
MAP a ramp can still outrun Newton, in which case the guess is restored and
the ramp is re-run with twice as many steps, twice: the ladder the
transient's cold start climbs (``simulation.py``, ``ISMIP7_CONTINUATION_STEPS``).
The inversion's warm start climbs the same ladder through this module.
"""
import numpy as np
from firedrake.exceptions import ConvergenceError


def ladder(base_steps):
    r"""The step counts tried in turn: ``base_steps``, then twice, then four
    times as many."""
    base = int(base_steps)
    if base < 1:
        raise ValueError(f"a continuation needs at least one step, not {base}")
    return (base, 2 * base, 4 * base)


def ramp_exponents(solve, z, n_flow, m_slide, n_target, m_target, steps_ladder,
                   report=print):
    r"""Ramp ``n_flow`` and ``m_slide`` from 1 to their targets, solving after
    each step, climbing ``steps_ladder`` on divergence.

    ``solve(attempt=, step=, steps=, t=)`` performs one solve for the current
    exponents and raises :class:`firedrake.exceptions.ConvergenceError` when
    it diverges; ``z`` is the state it solves for, restored to its entry value
    before each retry. Returns ``(attempt, steps)`` of the ramp that reached
    the targets; re-raises the last divergence when every rung fails, with
    ``z`` and the exponents back at their entry values.
    """
    z_init = z.copy(deepcopy=True)
    rungs = list(steps_ladder)
    last_error = None
    for attempt, steps in enumerate(rungs, 1):
        try:
            for step, t in enumerate(np.linspace(0.0, 1.0, steps), 1):
                n_flow.assign(1.0 + t * (n_target - 1.0))
                m_slide.assign(1.0 + t * (m_target - 1.0))
                solve(attempt=attempt, step=step, steps=steps, t=float(t))
            return attempt, steps
        except ConvergenceError as err:
            last_error = err
            z.assign(z_init)
            n_flow.assign(1.0)
            m_slide.assign(1.0)
            if attempt < len(rungs):
                report(f"  Continuation diverged at {steps} steps; "
                       f"restarting with {rungs[attempt]}...")
            else:
                report(f"  Continuation diverged at {steps} steps - giving up.")
    raise last_error
