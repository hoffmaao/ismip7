r"""Global (MPI-collective) statistics over Firedrake Functions.

``f.dat.data_ro`` is the OWNED slice on the calling rank. Any statistic taken
from it -- mean, min, max, sum, count -- is therefore rank-local, and under
``mpiexec`` it silently becomes a partition-dependent number. That is easy to
miss because it is correct on one rank, which is how most of this codebase was
developed and tested.

The defect was found in nine-plus places on this project. Most were printed
diagnostics (the theta/phi range after "Optimization finished", the C_w0 and
A_prior ranges, the dH/dt constraint's cell count, the MEaSUREs sigma
percentiles), but four were FUNCTIONAL:

- ``u_c = Constant(float(u_speed.dat.data_ro.mean()))`` in three drivers, which
  feeds the Weertman sliding coefficient. Each rank built friction from a
  different reference speed, so the same physical location got a different
  coefficient depending on which rank owned it.
- ``L = float(np.max(np.abs(coords.dat.data)))`` in the eikonal solver, which
  then divides the MESH COORDINATES by that factor. Each rank rescaled the
  geometry differently, so an element straddling a partition boundary was
  assembled from two different coordinate scalings.
- the sensitivity ``dJ_deps`` in gl_sensitivity.py, a rank-local dot product
  reported as the published Gt/yr-per-metre result.

**Two rules, and the second is not optional.**

1. Any statistic from ``.dat.data_ro`` is rank-local until reduced.
2. A collective is only safe where the surrounding control flow is itself
   rank-invariant. Wrapping a rank-local statistic in an allreduce inside a
   loop whose length differs per rank converts a wrong number into a DEADLOCK
   -- which happened here, in a ``sorted(..., key=lambda n: global_sum(...))``
   over a dict whose keys were built from rank-local masks, so ranks issued
   different numbers of collectives. Agree the iteration set across ranks
   first, then reduce in a deterministic order.
"""

import numpy as np


def _comm_of(f, comm):
    if comm is not None:
        return comm
    return f.function_space().mesh().comm


def global_sum(f, comm=None):
    r"""Sum over all owned dofs on every rank. Accepts a Function or an array
    (an array requires an explicit ``comm``, having no mesh to ask)."""
    data = f.dat.data_ro if hasattr(f, "dat") else np.asarray(f)
    if comm is None:
        if not hasattr(f, "dat"):
            raise TypeError(
                "global_sum on a bare array needs an explicit comm= "
                "(an ndarray carries no mesh to derive one from)"
            )
        comm = _comm_of(f, None)
    return float(comm.allreduce(float(np.sum(data))))


def global_dof_count(f, comm=None):
    r"""Number of owned SCALAR dofs across all ranks.

    ``.dat.data_ro.size`` on a vector-valued Function is dofs x components, so
    this divides by the block size -- passing a velocity to a naive ``.size``
    sum overcounts by the geometric dimension.
    """
    comm = _comm_of(f, comm)
    bs = getattr(f.function_space().dof_dset, "cdim", 1) or 1
    return int(comm.allreduce(int(f.dat.data_ro.size // bs)))


def global_mean(f, comm=None):
    r"""Area-agnostic mean over all owned dofs on every rank.

    This is the dof mean, not an area-weighted mean; use it where the old code
    used ``.dat.data_ro.mean()`` so behaviour matches on one rank.
    """
    comm = _comm_of(f, comm)
    data = f.dat.data_ro
    total = float(comm.allreduce(float(np.sum(data))))
    count = int(comm.allreduce(int(data.size)))
    return total / max(count, 1)


def global_range(f, comm=None):
    r"""``(min, max)`` over all owned dofs on every rank.

    Ranks owning no dofs contribute the identity element rather than a magic
    sentinel, so an empty partition cannot poison the result.
    """
    from mpi4py import MPI

    comm = _comm_of(f, comm)
    data = f.dat.data_ro if hasattr(f, "dat") else np.asarray(f)
    lo = float(np.min(data)) if data.size else np.inf
    hi = float(np.max(data)) if data.size else -np.inf
    return (float(comm.allreduce(lo, op=MPI.MIN)),
            float(comm.allreduce(hi, op=MPI.MAX)))


def global_max_abs(arr, comm):
    r"""``max(|arr|)`` across ranks, for raw arrays such as mesh coordinates.

    Takes an explicit comm: the caller is typically holding
    ``mesh.coordinates.dat.data``, which is not a Function.
    """
    from mpi4py import MPI

    a = np.asarray(arr)
    loc = float(np.max(np.abs(a))) if a.size else -np.inf
    return float(comm.allreduce(loc, op=MPI.MAX))


def global_count(mask, comm):
    r"""Number of True entries across ranks. Takes an explicit comm because a
    bare boolean ndarray has no mesh to derive one from."""
    m = np.asarray(mask)
    return int(comm.allreduce(int(np.count_nonzero(m))))
