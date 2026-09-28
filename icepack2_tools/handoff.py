r"""handoff.py - keeping one objective across an inversion's restarts.

A chained inversion warm-starts each link from the previous link's periodic
checkpoint. The optimiser only continues the same minimisation if the link
minimises the same objective from the point the previous link had accepted:

* every setting the objective depends on is recorded in the checkpoint and
  compared here on the warm start (:func:`objective_mismatches`); a link
  that would change one refuses under ``ISMIP7_WARM_START_STRICT=1``;
* the checkpoint holds the controls of the last ACCEPTED iterate, not of the
  last evaluated trial point of the line search
  (:func:`accepted_evaluation`);
* the checkpoint records the objective at that iterate, and the next link's
  first evaluation must reproduce it (:func:`handoff_gap`).

Rice (26 Sep 2026): "between restarts I think we need to check to make
sure the optimizer doesn't change."
"""
import math

# Checkpoint attributes the objective depends on. Auto-derived weights
# (log_vel_weight, prior_sigma_alpha) are FROZEN from the warm start rather
# than re-derived, so they appear here as plain values.
OBJECTIVE_KEYS = (
    "misfit_norm", "misfit_scale", "log_vel_weight", "log_vel_eps",
    "dhdt_weight", "dhdt_net_sigma",
    "prior_form", "gamma_theta", "gamma_phi",
    "prior_sigma_theta", "prior_sigma_phi", "prior_rho",
    "prior_sigma_alpha", "prior_rho_theta", "friction_c_ref",
    "friction_control", "friction", "n_flow", "geometry_space",
    "friction_anchor_length", "lake_ice_base", "fluidity_prior_origin",
    "grad_precond", "subelement_friction", "exact_front", "fluidity_control",
)

# Recorded with the controls: the objective at the checkpointed iterate.
OBJECTIVE_RECORD_KEYS = (
    "objective_total", "objective_misfit", "objective_reg_theta",
    "objective_reg_phi", "objective_iteration",
)


def _same(a, b, rel_tol):
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)
    if math.isnan(fa) and math.isnan(fb):
        return True
    return math.isclose(fa, fb, rel_tol=rel_tol, abs_tol=0.0)


def objective_mismatches(recorded, current, rel_tol=1e-9):
    r"""The keys of ``current`` whose recorded value differs, as
    ``"key: recorded -> current"`` strings. Keys the record lacks are not
    mismatches (an older checkpoint), keys ``current`` lacks are ignored."""
    out = []
    for key in OBJECTIVE_KEYS:
        if key in recorded and key in current:
            if not _same(recorded[key], current[key], rel_tol):
                out.append(f"{key}: {recorded[key]!r} -> {current[key]!r}")
    return out


def accepted_evaluation(ring, f_val, rel_tol=1e-10):
    r"""The entry of ``ring`` (an iterable of ``(value, *payload)`` for the
    most recent evaluations) whose value is ``f_val`` to ``rel_tol``: the
    evaluation the optimiser accepted. None if no entry matches, e.g. after
    a restart of the ring."""
    best, best_gap = None, math.inf
    for entry in ring:
        value = float(entry[0])
        gap = abs(value - f_val) / max(abs(f_val), 1e-300)
        if gap < best_gap:
            best, best_gap = entry, gap
    return best if best_gap <= rel_tol else None


def handoff_gap(recorded_total, first_total):
    r"""Relative difference between the objective a checkpoint recorded and
    the value the warm-started link evaluates first: zero when the link
    continues the same minimisation from the same point."""
    return abs(float(first_total) - float(recorded_total)) / max(
        abs(float(recorded_total)), 1e-300)
