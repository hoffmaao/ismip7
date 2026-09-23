r"""Reject a diagnostic solve that converged onto a speed no ice can have.

The step loop asks the solver whether it converged and nothing else, so a
solve that converges onto a runaway is accepted and carried forward. That is
how the Lambert/Amery trough destroys a run: SNES reason=2, function norm
below 1, peak speed 2.5e6 m/yr, step accepted. The caller raises a violation
as a convergence failure, which hands the step to the rescue ladder and then
the subcycles.
"""

from firedrake import dot, sqrt

from icepack2_tools.mpi_stats import global_extreme_location

# A node the soft limiter pins settles where its extra drag
# k_lim * (|u| - u_lim) balances the driving stress, so necessarily above
# u_lim: 2.01e4 to 2.1e4 m/yr at the default k_lim and u_lim = 2e4. A bound
# equal to u_lim would veto exactly the rungs whose job is to pin a runaway.
LIMITER_BOUND_FACTOR = 3.0


def limiter_speed_bound(bound, u_lim):
    r"""The bound for a rescue rung that runs the soft speed limiter.

    The ordinary path keeps the configured ``bound`` so it still catches the
    trough the step it ignites. A limiter rung takes the larger of ``bound``
    and ``LIMITER_BOUND_FACTOR * u_lim``, which sits above where a pinned node
    settles and far below the runaway. ``0`` stays disabled.
    """
    if bound <= 0.0:
        return 0.0
    return max(bound, LIMITER_BOUND_FACTOR * u_lim)


def speed_bound_violation(u, speed, coordinates, bound):
    r"""``(u_max, (x, y))`` when the peak of ``|u|`` exceeds ``bound``, else
    ``None``.

    ``speed`` is a scalar work Function the speed is interpolated into and
    ``coordinates`` holds its dof coordinates. Collective over the mesh
    communicator. ``bound <= 0`` disables the check before any interpolation
    or reduction, which the step loop would otherwise pay after every solve.
    """
    if bound <= 0.0:
        return None
    speed.interpolate(sqrt(dot(u, u)))
    u_max, u_max_xy = global_extreme_location(speed, coordinates, mode="max")
    if u_max <= bound:
        return None
    return u_max, u_max_xy
