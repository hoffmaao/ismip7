#!/bin/bash
# ISMIP7 icepack2 inversion on Rice NOTS, with a self-resuming chain.
#
# Submit through nots_submit_inversions.sh, which sets the mesh, sampling and
# output name; or directly:
#   sbatch -p long -C sapphirerapids --mem=240G -N1 -n32 --time=3-00:00:00 \
#     --export=ALL,ISMIP7_LC=2000,ISMIP7_LC_COARSE=5000,ISMIP7_BUFFER_M=0,\
#ISMIP7_MESH=...,ISMIP7_RASTER_SAMPLE=vertex,ISMIP7_MAP_OUT=... \
#     antarctica/scripts/batch_runners/nots_inversion.sbatch
#
# CHAIN. A 2 km inversion can outlast a wall limit. The inversion checkpoints
# ISMIP7_MAP_OUT every 20 iterates and can warm-start from it, so this script
# queues its own successor FIRST (--dependency=afterany, so it runs whether
# this job converged, timed out or lost its node), and each link exits at once
# if the "<map>.done" marker exists. Depth is capped by ISMIP7_CHAIN_MAX. The
# successor copies this job's partition, constraint, memory, time and task
# layout from scontrol, so a chain never changes rank count mid-run -- the
# inversion now refuses a warm start whose dof ordering differs anyway.
#
# Defaults write inversion_icepack2_rc_n3_dg0_logvelnet_2500.h5 under the
# settings the 2500 m result came from: the sigma-normalised velocity chi^2
# plus ISSM's logarithmic misfit, the pointwise dH/dt term, and the integrated
# net mass-balance term (off by default in the repo) that takes grounding-line
# discharge from 1628 to 2211 Gt/yr. The friction law is the one deliberate
# difference: nots_env.sh defaults ISMIP7_FRICTION to regularized_coulomb.

#SBATCH -J ismip7_inv
#SBATCH -p long
#SBATCH -o logs/ismip7_inv_%j.out
#SBATCH -e logs/ismip7_inv_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=END,FAIL

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p logs
. antarctica/scripts/batch_runners/nots_env.sh
nots_activate
nots_banner

export ISMIP7_MISFIT_NORM="${ISMIP7_MISFIT_NORM:-sigma}"
export ISMIP7_SIGMA_U_FLOOR="${ISMIP7_SIGMA_U_FLOOR:-3}"
export ISMIP7_LOG_VEL_WEIGHT="${ISMIP7_LOG_VEL_WEIGHT:-auto}"
export ISMIP7_DHDT_WEIGHT="${ISMIP7_DHDT_WEIGHT:-1}"
export ISMIP7_DHDT_NET_SIGMA="${ISMIP7_DHDT_NET_SIGMA:-10}"
export ISMIP7_GAMMA_THETA="${ISMIP7_GAMMA_THETA:-1e5}"
export ISMIP7_GAMMA_PHI="${ISMIP7_GAMMA_PHI:-1e5}"
export ISMIP7_MAXITER="${ISMIP7_MAXITER:-200}"
export ISMIP7_RASTER_SAMPLE="${ISMIP7_RASTER_SAMPLE:-vertex}"
# Name the output explicitly: without it every run, including a smoke test,
# writes the production filename and can silently replace a converged map.
# ISMIP7_MAP_DEFAULT is named from the law, so ISMIP7_FRICTION=budd writes a
# _budd file rather than resuming from and then overwriting the RC one. It is
# the same path nots_projection.sbatch loads.
export ISMIP7_MAP_OUT="${ISMIP7_MAP_OUT:-$ISMIP7_MAP_DEFAULT}"
DONE="${ISMIP7_MAP_OUT}.done"
LOG="logs/ismip7_inv_${SLURM_JOB_ID}.out"

echo "    map out  $ISMIP7_MAP_OUT"
echo "    sampling $ISMIP7_RASTER_SAMPLE  net sigma $ISMIP7_DHDT_NET_SIGMA  log-vel $ISMIP7_LOG_VEL_WEIGHT  maxiter $ISMIP7_MAXITER"

if [ -f "$DONE" ]; then
    echo "chain: $DONE exists, inversion already finished; nothing to do"
    exit 0
fi

# --- queue the successor before doing anything slow -----------------------
DEPTH="${ISMIP7_CHAIN_DEPTH:-0}"
MAXD="${ISMIP7_CHAIN_MAX:-4}"
if [ "$DEPTH" -lt "$MAXD" ]; then
    INFO=$(scontrol show job "$SLURM_JOB_ID")
    feat=$(echo "$INFO" | grep -oE 'Features=[^ ]+' | cut -d= -f2); [ "$feat" = "(null)" ] && feat=""
    tlim=$(echo "$INFO" | grep -oE 'TimeLimit=[^ ]+' | cut -d= -f2)
    memn=$(echo "$INFO" | grep -oE 'MinMemoryNode=[^ ]+' | cut -d= -f2)
    NEXT=$(sbatch --parsable --dependency=afterany:"$SLURM_JOB_ID" \
        -p "$SLURM_JOB_PARTITION" ${feat:+-C "$feat"} --time="$tlim" --mem="$memn" \
        -N "$SLURM_JOB_NUM_NODES" -n "$SLURM_NTASKS" --cpus-per-task="${SLURM_CPUS_PER_TASK:-1}" \
        --export=ALL,ISMIP7_CHAIN_DEPTH=$((DEPTH + 1)) \
        antarctica/scripts/batch_runners/nots_inversion.sbatch) \
      && echo "chain: successor $NEXT queued (depth $((DEPTH + 1)) of $MAXD)" \
      || echo "chain: WARNING could not queue a successor"
else
    echo "chain: depth $DEPTH reached ISMIP7_CHAIN_MAX=$MAXD, no successor"
fi

# --- resume if a checkpoint exists ------------------------------------------
if [ -f "$ISMIP7_MAP_OUT" ]; then
    export ISMIP7_WARM_START="$ISMIP7_MAP_OUT"
    echo "chain: warm-starting theta/phi from the periodic checkpoint $ISMIP7_MAP_OUT"
fi

srun -n "$SLURM_NTASKS" python -u antarctica/scripts/inversion_icepack2.py
rc=$?
echo "=== $(date -Is)  inversion exit ${rc} ==="

# Finished for good only if L-BFGS-B returned for a reason other than its
# iteration cap; an ITERATIONS REACHED LIMIT return is what a successor is for.
if grep -q 'Optimization finished' "$LOG" && ! grep -q 'ITERATIONS REACHED LIMIT' "$LOG"; then
    touch "$DONE"
    echo "chain: finished ($(grep -o 'Optimization finished.*' "$LOG" | tail -1)); marker $DONE"
fi
exit $rc
