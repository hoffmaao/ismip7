#!/bin/bash
# Prove the NOTS Firedrake build works ACROSS RANKS under srun, and that the
# parallel partitioners actually partition. The build job's own check ran
# `srun -n 4` inside a 1-task allocation and so failed for every type,
# `simple` included, which says nothing about the build. This job has 4 tasks.
# Resources (partition, tasks, memory, wall time, account, constraint) come
# from submit.sh and this cluster's sites/<name>.sh, never from directives
# here: a directive cannot read a site file, and a header that disagrees
# with the command line is rejected outright by Slurm.
#SBATCH -o logs/fd_verify_%j.out
#SBATCH -e logs/fd_verify_%j.err
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
. antarctica/scripts/batch_runners/site_env.sh
ismip7_activate
ismip7_banner
export PYOP2_CACHE_DIR="${PYOP2_CACHE_DIR:-$HOME/.pyop2_cache_verify}"
echo "firedrake: $(python -c 'import firedrake,sys; print(firedrake.__version__)' 2>/dev/null | tail -1)"
cat > "$PYOP2_CACHE_DIR/_part.py" <<'PYEOF'
from firedrake import *
from firedrake.petsc import PETSc
import sys
PETSc.Options()["petscpartitioner_type"] = sys.argv[1]
m = UnitSquareMesh(80, 80, distribution_parameters={"partition": True})
V = FunctionSpace(m, "CG", 1)
owned, total = V.dof_dset.size, V.dof_dset.total_size
g = m.comm.allreduce(total - owned) / max(m.comm.allreduce(owned), 1)
if m.comm.rank == 0:
    print(f"    {sys.argv[1]:9s} ranks={m.comm.size} ghost/owned = {g:.4f}")
PYEOF
for p in parmetis ptscotch simple; do
    srun -n 4 python "$PYOP2_CACHE_DIR/_part.py" "$p" 2>"$PYOP2_CACHE_DIR/_err_$p.txt" | grep ghost \
        || { echo "    $p FAILED:"; grep -iE 'error|Traceback' "$PYOP2_CACHE_DIR/_err_$p.txt" | head -3; }
done
echo; echo "=== real mesh: antarctica_64000_2500.msh at 4 ranks, ptscotch vs simple ==="
for p in ptscotch simple; do
    srun -n 4 python -u antarctica/scripts/partition_probe.py antarctica/mesh/antarctica_64000_2500.msh \
        -petscpartitioner_type "$p" 2>/dev/null | grep -E 'ranks=' | sed "s/^/    $p: /"
done
echo "=== $(date -Is) verify done ==="
