r"""Global (cross-rank) statistics of a distributed Function.

Anything read off ``Function.dat.data_ro`` covers only the dofs the calling
rank OWNS. ``PETSc.Sys.Print`` emits from rank 0 alone, so printing such a
statistic reports rank 0's slice of the field: the number reads narrower than
the field really is and changes with the rank count. Worse, using one as a
scalar in a form (a normalization constant, say) gives each rank a different
coefficient, so the assembled residual and Jacobian disagree across ranks and
the answer depends on the partition.

Reduce first, through these helpers. Every one of them is collective: call it
on all ranks or not at all.

See AGENTS.md section 3.
"""

import numpy as np
from mpi4py import MPI


def _comm_of(f, comm):
    return f.function_space().mesh().comm if comm is None else comm


def global_range(f, comm=None):
    r"""``(min, max)`` of a Function over all ranks."""
    d = np.asarray(f.dat.data_ro)
    comm = _comm_of(f, comm)
    return (comm.allreduce(float(d.min()) if d.size else np.inf, op=MPI.MIN),
            comm.allreduce(float(d.max()) if d.size else -np.inf, op=MPI.MAX))


def global_max(f, comm=None):
    r"""Maximum of a Function over all ranks."""
    d = np.asarray(f.dat.data_ro)
    return _comm_of(f, comm).allreduce(
        float(d.max()) if d.size else -np.inf, op=MPI.MAX)


def global_mean(f, comm=None):
    r"""Mean of a Function over the owned dofs of all ranks."""
    d = np.asarray(f.dat.data_ro)
    comm = _comm_of(f, comm)
    total = comm.allreduce(float(d.sum()), op=MPI.SUM)
    n = comm.allreduce(int(d.size), op=MPI.SUM)
    return total / max(n, 1)


def global_size(f, comm=None):
    r"""Total number of owned dofs across all ranks."""
    return _comm_of(f, comm).allreduce(int(f.dat.data_ro.size), op=MPI.SUM)


def global_count(mask, comm):
    r"""Number of True entries of a rank-local boolean array, over all ranks."""
    return comm.allreduce(int(np.count_nonzero(mask)), op=MPI.SUM)
