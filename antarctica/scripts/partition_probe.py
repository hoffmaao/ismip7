"""Measure mesh partition quality vs rank count.

Wall-clock scaling cannot be measured on a machine at load 147/80, but
partition quality is a property of the graph and the partitioner alone, so it
is immune to whatever else is running. It is also the thing most likely to
decide how this code scales on a cluster: with no parallel graph partitioner
in the PETSc build, Firedrake falls back to PETSc's "simple" partitioner,
which cuts the DMPlex by cell index rather than by locality.

Reports, for a CG1 field (the control/velocity space):
  imbalance     max owned cells / mean owned cells   (1.0 is perfect)
  ghost ratio   ghost dofs / owned dofs, per rank    (communication volume)
A locality-preserving partition of a 2D mesh has ghost/owned ~ C/sqrt(n_local),
so it should FALL as ranks rise for fixed mesh size only slowly, and stay
well under 1. A ratio near or above 1 means most of each rank's data is halo.
"""
import sys
import numpy as np
from firedrake import Mesh, FunctionSpace
from firedrake.petsc import PETSc
from mpi4py import MPI

fn = sys.argv[1]
comm = MPI.COMM_WORLD
mesh = Mesh(fn, comm=comm)
Q = FunctionSpace(mesh, "CG", 1)
DG = FunctionSpace(mesh, "DG", 0)

owned_dof = Q.dof_dset.size
total_dof = Q.dof_dset.total_size
ghost_dof = total_dof - owned_dof
owned_cells = DG.dof_dset.size

oc = np.array(comm.allgather(owned_cells), dtype=float)
od = np.array(comm.allgather(owned_dof), dtype=float)
gd = np.array(comm.allgather(ghost_dof), dtype=float)

if comm.rank == 0:
    n = comm.size
    ratio = gd / np.maximum(od, 1)
    print(f"  ranks={n:3d}  cells={int(oc.sum()):>9,}  dofs={int(od.sum()):>9,}"
          f"  imbalance={oc.max()/oc.mean():5.2f}"
          f"  ghost/owned mean={ratio.mean():6.3f} max={ratio.max():6.3f}"
          f"  halo dofs={int(gd.sum()):>9,} ({gd.sum()/od.sum()*100:5.1f}% of owned)")
    sys.stdout.flush()
