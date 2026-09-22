r"""Cross-mesh transfer of checkpoint fields with an explicit fill.

A forward that names a MAP with ``ISMIP7_INVERSION`` and a compute mesh with
``ISMIP7_MESH`` interpolates the MAP's continuous fields onto the compute
mesh. Firedrake's cross-mesh ``interpolate`` can only evaluate the source
where the source mesh exists; a target dof outside it is a "missing dof".
The stock behaviour writes ``0.0`` there, which is harmless for the log
controls (theta = phi = 0 is the prior) and fatal for the fluidity prior: the
2 km MAPs were inverted on a buffer-0 mesh and the production mesh carries a
20 km ocean buffer, so every dof in that ring is missing, and
``A_eff = A_prior * exp(phi) = 0`` there zeroes the dislocation term, the
``lin_reg`` regularizer and the ``alpha_gl`` collar in
``dual_friction.build_rc_residual`` at once: a singular membrane block.

``interpolate_with_fill`` makes the fill a stated choice and counts it.
"""
import numpy as np
from firedrake import Function

from .mpi_stats import global_count, global_size


def interpolate_with_fill(target, source, fill, comm=None):
    r"""Interpolate ``source`` into ``target`` across meshes; dofs of ``target``
    outside the source mesh take ``fill``.

    ``fill`` is a float, or a Function on ``target``'s space whose values are
    taken where the source has none (the raster-sampled velocity_obs, say).
    Returns ``(n_missing, n_total)``, both reduced over ranks and counting
    dofs once (owned dofs only, whatever the value shape). A same-mesh call
    is a plain interpolate and reports no missing dofs.
    """
    comm = comm if comm is not None else target.comm
    total = global_size(target, comm)
    if source.function_space().mesh() is target.function_space().mesh():
        target.interpolate(source)
        return 0, total
    target.interpolate(
        source, allow_missing_dofs=True, default_missing_val=np.nan
    )
    data = target.dat.data
    missing = np.isnan(data.reshape(data.shape[0], -1)).any(axis=1)
    if isinstance(fill, Function):
        if fill.function_space() != target.function_space():
            raise ValueError(
                "the fill Function must live on the target's function space"
            )
        data[missing] = fill.dat.data_ro[missing]
    else:
        data[missing] = float(fill)
    n_missing = global_count(missing, comm)
    if global_count(np.isnan(data.reshape(data.shape[0], -1)).any(axis=1), comm):
        raise RuntimeError(
            f"{target.name()}: NaN left after filling missing dofs; the fill "
            "itself carries NaN"
        )
    return n_missing, total
