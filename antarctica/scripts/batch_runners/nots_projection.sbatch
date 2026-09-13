#!/bin/bash
# ISMIP7 forward run (control, historical or projection) on Rice NOTS deepsC.
# SELF-CHAINING: a 285-year projection is about five days at 2500 m and will
# not fit any ordinary queue, so this script resubmits itself until the run
# reaches ISMIP7_T_END, resuming from its own newest checkpoint each time.
#
#   sbatch --export=ALL,ISMIP7_EXPERIMENT=control \
#          antarctica/scripts/batch_runners/nots_projection.sbatch
#
# ISMIP7_EXPERIMENT selects the driver:
#   control            antarctica/scripts/control/run.py
#   ssp126_cesm_waccm  antarctica/scripts/projections/ssp126_cesm_waccm.py
#   ocx                antarctica/scripts/projections/ocx.py
#   hist_cesm_waccm    antarctica/scripts/historical/cesm_waccm.py
#   hist_mri_esm2      antarctica/scripts/historical/mri_esm2.py
#
# Sizing, measured on Cascade Lake at the same clock: 155 s per 0.1 yr step at
# 2500 m on 12 ranks, i.e. ~26 min per simulated year, ~5 core-hours per year.
# 24 h of wall time therefore buys roughly 55 simulated years, so a 285-year
# projection is about six chained jobs.
#
# NOTE (measured 2026-09-09): deepsC has NO wall-time limit, so on deepsC the
# 24 h below is a self-imposed courtesy on a shared condo and the chaining is
# optional. The chaining is what makes this runnable on `commons` (1-day limit)
# and `long` (3-day limit), and it is also the safety net that turns a node
# failure into a resume rather than a lost week. Keep it on.

#SBATCH -J ismip7_fwd
#SBATCH -p commons
#SBATCH -C cascadelake
#SBATCH -o logs/ismip7_fwd_%j.out
#SBATCH -e logs/ismip7_fwd_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p logs
. antarctica/scripts/batch_runners/nots_env.sh
nots_activate
nots_banner

EXPERIMENT="${ISMIP7_EXPERIMENT:-control}"
case "$EXPERIMENT" in
    control)           DRIVER=antarctica/scripts/control/run.py ;;
    ssp126_cesm_waccm) DRIVER=antarctica/scripts/projections/ssp126_cesm_waccm.py ;;
    ocx)               DRIVER=antarctica/scripts/projections/ocx.py ;;
    hist_cesm_waccm)   DRIVER=antarctica/scripts/historical/cesm_waccm.py ;;
    hist_mri_esm2)     DRIVER=antarctica/scripts/historical/mri_esm2.py ;;
    *) echo "ERROR: unknown ISMIP7_EXPERIMENT='$EXPERIMENT'" >&2; exit 2 ;;
esac

# --- forward configuration ---------------------------------------------
export ISMIP7_DT="${ISMIP7_DT:-0.1}"
export ISMIP7_T_END="${ISMIP7_T_END:-2300}"
export ISMIP7_CHECKPOINT_EVERY_YR="${ISMIP7_CHECKPOINT_EVERY_YR:-5}"
export ISMIP7_KEEP_CHECKPOINTS="${ISMIP7_KEEP_CHECKPOINTS:-3}"
# The control configuration the protocol asks for: calving held at end-of-2014
# conditions, with the apparent mass-balance reference. Projections that should
# calve freely must override ISMIP7_CALVING and unset ISMIP7_FIXED_FRONT.
export ISMIP7_APPARENT_MB="${ISMIP7_APPARENT_MB:-1}"
export ISMIP7_FIXED_FRONT="${ISMIP7_FIXED_FRONT:-1}"
# Resume from this experiment's own newest checkpoint. This is what makes the
# chain below safe to run unattended.
export ISMIP7_AUTO_RESUME="${ISMIP7_AUTO_RESUME:-1}"

echo "    experiment $EXPERIMENT -> $DRIVER"
echo "    dt=$ISMIP7_DT t_end=$ISMIP7_T_END auto_resume=$ISMIP7_AUTO_RESUME" \
     "apparent_mb=$ISMIP7_APPARENT_MB fixed_front=$ISMIP7_FIXED_FRONT"

# --- run, leaving time to resubmit cleanly ------------------------------
# Stop the model a little before the wall clock so the final checkpoint is
# written and the next job in the chain starts from a complete state.
srun -n "${SLURM_NTASKS}" python -u "$DRIVER"
rc=$?
echo "=== $(date -Is)  forward exit ${rc} ==="

# --- chain --------------------------------------------------------------
# Resubmit only when the driver exited cleanly AND the run has not reached
# t_end. A non-zero exit is a real failure and must not be retried blindly:
# the July grounding-line blow-up looked exactly like a run that "just needed
# more time", and chaining through it would have burned days.
CHAIN="${ISMIP7_CHAIN:-1}"
if [ "$rc" -ne 0 ]; then
    echo "chain: stopping, driver exited ${rc}"
    exit $rc
fi
if [ "$CHAIN" != "1" ]; then
    echo "chain: disabled (ISMIP7_CHAIN=$CHAIN)"
    exit 0
fi

# The drivers print the reached year; ask the timeseries instead, which is the
# state the next job would resume from.
LAST_YEAR=$(python - <<'PY'
import csv, glob, os, sys
tag = os.environ.get("ISMIP7_RUN_TAG", "")
pats = sorted(glob.glob("antarctica/results/*timeseries.csv"), key=os.path.getmtime)
if not pats:
    sys.exit(0)
rows = list(csv.DictReader(open(pats[-1])))
print(rows[-1]["year"] if rows else "")
PY
)
T_END="${ISMIP7_T_END}"
if [ -z "$LAST_YEAR" ]; then
    echo "chain: no timeseries found, not resubmitting"
    exit 0
fi
done_yet=$(python -c "print(1 if float('$LAST_YEAR') >= float('$T_END') - 1e-6 else 0)")
if [ "$done_yet" = "1" ]; then
    echo "chain: reached t_end ${T_END} at ${LAST_YEAR}, finished"
    exit 0
fi

echo "chain: at ${LAST_YEAR} of ${T_END}, resubmitting"
sbatch --dependency=afterok:"${SLURM_JOB_ID}" --export=ALL \
       antarctica/scripts/batch_runners/nots_projection.sbatch
