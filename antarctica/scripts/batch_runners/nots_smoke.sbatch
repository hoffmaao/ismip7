#!/bin/bash
# Smoke test: the 32 km inversion for 2 iterates on 4 ranks. Exercises the
# whole stack on NOTS (firedrake, icepack, icepack2, tlm_adjoint, rasterio,
# every data file) in a few minutes, so a multi-day job cannot die on a
# missing file or module hours after it queued.
#   sbatch antarctica/scripts/batch_runners/nots_smoke.sbatch
#SBATCH -J ismip7_smoke
#SBATCH -p scavenge
#SBATCH -C cascadelake
#SBATCH -o logs/ismip7_smoke_%j.out
#SBATCH -e logs/ismip7_smoke_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=24G
#SBATCH --time=00:45:00
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
# Defaults are the 32 km smoke; every value can be overridden through
# sbatch --export=ALL,... to smoke-test another mesh or sampling, e.g. two
# iterates of the 2 km cell-mean case on 32 Sapphire Rapids ranks. They MUST
# be set before nots_env.sh is sourced: it fills ISMIP7_LC / ISMIP7_MESH with
# the 2500 m production defaults, after which a ${VAR:-...} fallback here never
# fires (job 1300996 silently ran the 2500 m mesh on 4 ranks that way).
export ISMIP7_LC="${ISMIP7_LC:-32000}" ISMIP7_LC_COARSE="${ISMIP7_LC_COARSE:-320000}" ISMIP7_BUFFER_M="${ISMIP7_BUFFER_M:-0}"
export ISMIP7_MESH="${ISMIP7_MESH:-/projects/ah301/ismip7/antarctica/mesh/antarctica_320000_32000.msh}"
. antarctica/scripts/batch_runners/nots_env.sh
nots_activate
nots_banner
export ISMIP7_MISFIT_NORM=sigma ISMIP7_SIGMA_U_FLOOR=3 ISMIP7_LOG_VEL_WEIGHT=auto
export ISMIP7_DHDT_WEIGHT=1 ISMIP7_DHDT_NET_SIGMA=10 ISMIP7_GAMMA_THETA=1e5 ISMIP7_GAMMA_PHI=1e5
export ISMIP7_MAXITER="${ISMIP7_MAXITER:-2}" ISMIP7_RASTER_SAMPLE="${ISMIP7_RASTER_SAMPLE:-vertex}"
export ISMIP7_MAP_OUT="${ISMIP7_MAP_OUT:-$ISMIP7_REPO/antarctica/results/smoke_${ISMIP7_LC}_${ISMIP7_LC_COARSE}_${ISMIP7_RASTER_SAMPLE}.h5}"
mkdir -p "$(dirname "$ISMIP7_MAP_OUT")"
srun -n "$SLURM_NTASKS" python -u antarctica/scripts/inversion_icepack2.py
rc=$?
echo "=== $(date -Is) smoke exit $rc ==="
exit $rc
