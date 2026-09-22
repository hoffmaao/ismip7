#!/bin/bash
# Partition-quality probe. RUN THIS FIRST, before spending queue time on a
# timing or production job.
#
#   antarctica/scripts/batch_runners/submit.sh probe
#
# It goes to this cluster's ISMIP7_PART_DEBUG queue (at Rice that is
# `scavenge`: 1-hour limit, preemptible). The probe is short and idempotent, so
# preemption costs nothing and a debug queue is far shorter than a production
# one.
#
# It only distributes the mesh: no solve, no data, about a minute per rank
# count, and immune to whatever else is on the node. What it reports is the
# ghost-to-owned dof ratio per rank, which is the one number that says whether
# the Firedrake/PETSc build on this machine can be parallelised at all.
#
# WHAT TO EXPECT. A locality-preserving partition of a 2D mesh gives
# ghost/owned of roughly 4*sqrt(P/N), so a few hundredths here, falling slowly
# with rank count. On a build WITHOUT ParMETIS or PTScotch, PETSc falls back to
# its "simple" partitioner, which cuts the mesh by point index rather than by
# locality; measured on the workstation build that gives
#
#     ranks    2      4      8     16     32
#     mean  2.97   16.0   59.6  187.1  287.5      ghost/owned per rank
#     halo   92%   182%   313%   480%   653%      of owned dofs
#
# i.e. 284x the ideal at 2 ranks and 6330x at 16, growing as about P^1.6. If
# this job reports anything in that range, the build is missing a parallel
# partitioner and every scaling number taken on it would measure the
# partitioner rather than the model. Rebuild PETSc with --download-parmetis
# --download-ptscotch before going further.

# Resources (partition, tasks, memory, wall time, account, constraint) come
# from submit.sh and this cluster's sites/<name>.sh, never from directives
# here: a directive cannot read a site file, and a header that disagrees
# with the command line is rejected outright by Slurm.
#SBATCH -o logs/ismip7_part_%j.out
#SBATCH -e logs/ismip7_part_%j.err

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p logs
. antarctica/scripts/batch_runners/site_env.sh
ismip7_activate
ismip7_banner

MESH="${ISMIP7_MESH}"
echo "    probing $MESH"
echo

# Only up to what this allocation actually holds: a step asking for more ranks
# than the site's ISMIP7_TASKS cannot be created, so it yields no measurement.
for n in 1 2 4 8 16 32; do
    [ "$n" -gt "${SLURM_NTASKS:-1}" ] && break
    ismip7_mpirun "$n" \
        python -u antarctica/scripts/partition_probe.py "$MESH" 2>&1 \
        | grep -iE "ranks=|error|traceback|signal|killed|fault"
done

echo
echo "=== $(date -Is)  probe done ==="
echo "Read the ghost/owned column. Hundredths: good. Tens or hundreds: the"
echo "build has no parallel graph partitioner; fix that before timing anything."
