r"""Cross-mesh transfer of a checkpoint field onto another mesh.

``simulation.setup_model`` interpolates a MAP's controls onto a different
compute mesh when ``ISMIP7_MESH`` names one (the timing matrix, a
production mesh finer than the inversion's). Firedrake locates each target
dof in the source mesh with the source mesh's ``tolerance``, 0.5 of the
reference cell by default, so a dof up to half a source cell OUTSIDE the
source domain is still assigned a cell and the field is EXTRAPOLATED there.
Measured (22 September 2026): the 2 km buffered0 Budd MAP transferred onto
``antarctica_50000_5000_buffered20000`` gave a fluidity prior in
[-164.7, 921.1] from a source range of [1.0, 783.7], and the condensed
diagnostic solve then aborted on a singular local block ("Getri throws
nonzero info"). The ocean buffer of the target lies outside the source
outline, which is exactly where the extrapolation happens.

:func:`strict_transfer` locates strictly, fills the dofs outside the
source with a value the caller chooses (the prior's floor for a fluidity,
zero for a log-deviation or an observation), and clamps the result to the
source's range, so a transferred field can never leave the range the
inversion produced.
"""
import numpy as np
from firedrake import Function, FunctionSpace

from .mpi_stats import global_count, global_range, global_size

#: Relative point-location tolerance for a strict transfer. Zero would let
#: floating-point error push a dof on the shared outline outside; this keeps
#: the outline and rejects anything a source cell does not contain.
STRICT_TOLERANCE = 1e-8


def outside_source(source_mesh, target_space):
    r"""Boolean per owned dof of ``target_space``: not inside any source cell.

    An indicator that is one on the whole source mesh interpolates to one at
    every located dof and to the missing-dof default, zero, elsewhere."""
    source_mesh.tolerance = STRICT_TOLERANCE
    one = Function(FunctionSpace(source_mesh, "CG", 1)).assign(1.0)
    el = target_space.ufl_element()
    scalar = FunctionSpace(target_space.mesh(), el.family(), el.degree())
    hit = Function(scalar).interpolate(
        one, allow_missing_dofs=True, default_missing_val=0.0)
    return hit.dat.data_ro < 0.5


def strict_transfer(source_field, target_space, fill=0.0, outside=None):
    r"""Interpolate ``source_field`` onto ``target_space`` without leaving
    the source domain or the source range.

    Returns ``(field, info)`` with ``info`` carrying ``n_outside`` and
    ``n_all`` (global dof counts), ``fill``, and the source and target
    ranges. A ``fill`` outside the source range is clipped into it, and
    ``info["fill"]`` is the value actually written. ``outside`` may be
    passed from :func:`outside_source` when several fields share one target
    space."""
    source_mesh = source_field.function_space().mesh()
    source_mesh.tolerance = STRICT_TOLERANCE
    field = Function(target_space, name=source_field.name())
    field.interpolate(source_field, allow_missing_dofs=True,
                      default_missing_val=0.0)
    if outside is None:
        outside = outside_source(source_mesh, target_space)
    n_out = global_count(outside, target_space.mesh().comm)
    n_all = global_size(field)
    src_lo, src_hi = global_range(source_field)
    fill = float(np.clip(fill, src_lo, src_hi))
    if n_out:
        field.dat.data[outside] = fill
    np.clip(field.dat.data, src_lo, src_hi, out=field.dat.data)
    tgt_lo, tgt_hi = global_range(field)
    info = {"n_outside": n_out, "n_all": n_all, "fill": fill,
            "source_range": (src_lo, src_hi), "target_range": (tgt_lo, tgt_hi)}
    return field, info


def describe(name, info):
    lo, hi = info["target_range"]
    slo, shi = info["source_range"]
    return (f"Transfer {name}: {info['n_outside']}/{info['n_all']} target dofs "
            f"outside the source domain -> {info['fill']:.3g}; target "
            f"[{lo:.4g}, {hi:.4g}] within source [{slo:.4g}, {shi:.4g}]")
