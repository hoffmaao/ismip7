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

Measured 2026-09-05 on antarctica_64000_2500.msh (576,801 cells, 293,006 CG1
dofs) with PETSc 3.25 built --download-metis only, i.e. with no ParMETIS,
PTScotch or Chaco, so Firedrake fell back to the PETSc simple partitioner.
Cell imbalance was 1.00 at every rank count, so what follows is communication
cost, not load imbalance:

  ranks   halo (% of owned dofs)   ghost/owned mean (worst rank)
      2                      92%                    2.97  (5.8)
      4                     182%                   16.0   (59)
      8                     313%                   59.6   (377)
     16                     480%                  187.1   (2343)
     32                     653%                  287.5   (4703)

A compact 2D partition of this mesh would give ghost/owned of roughly
4*sqrt(P/N), i.e. 0.01 to 0.03 over this range, so the measured halo is 284x
that at 2 ranks and 6330x at 16, and it grows roughly as P^1.6. Conclusion: no
scaling curve taken from this build measures the model, only the partitioner,
so rebuild PETSc with a parallel graph partitioner before reading one. Raw
output: antarctica/results/logs/partition_2500.log (results/ is gitignored).
"""
import sys
import numpy as np
from firedrake import Mesh, FunctionSpace
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
