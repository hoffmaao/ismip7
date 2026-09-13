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
# Point the forward at the MAP the inversion chain wrote. Without this it would
# derive the name through icepack2_tools.naming.map_basename, which has no
# `logvelnet` element, and abort on a missing file after the queue wait.
export ISMIP7_INVERSION="${ISMIP7_INVERSION:-$ISMIP7_MAP_DEFAULT}"
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
# An integer flag: 0 turns it off. --export=ALL,VAR=... cannot unset a
# variable, so 0 is the only way to ask a resubmission for a clean run.
# Normalised here so the banner below reports what the drivers will do.
case "${ISMIP7_AUTO_RESUME:-1}" in
    ''|0|00) export ISMIP7_AUTO_RESUME=0 ;;
    *)       export ISMIP7_AUTO_RESUME=1 ;;
esac

echo "    experiment $EXPERIMENT -> $DRIVER"
echo "    map in  $ISMIP7_INVERSION"
[ -f "$ISMIP7_INVERSION" ] \
    || echo "    WARNING: that MAP does not exist; run nots_inversion.sbatch first" \
            "or set ISMIP7_INVERSION"
echo "    dt=$ISMIP7_DT t_end=$ISMIP7_T_END auto_resume=$ISMIP7_AUTO_RESUME" \
     "apparent_mb=$ISMIP7_APPARENT_MB fixed_front=$ISMIP7_FIXED_FRONT"

# --- run, leaving time to resubmit cleanly ------------------------------
# Stop the model a little before the wall clock so the final checkpoint is
# written and the next job in the chain starts from a complete state. Slurm
# gives no warning at the limit, it just kills the step, so the budget has to
# be handed to the driver: ISMIP7_WALL_STOP_MIN is minutes of run time.
#
# Derived from THIS job's own TimeLimit, and deliberately not exported: the
# resubmit below uses --export=ALL, so an exported value would pin every later
# link in the chain to the first job's limit even when it lands on a partition
# with a shorter one.
wall_stop_min="${ISMIP7_WALL_STOP_MIN:-}"
tl=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
     | grep -oE 'TimeLimit=[^ ]+' | cut -d= -f2)
if [ -z "$wall_stop_min" ]; then
    wall_stop_min=$(python - "$tl" <<'PY'
import sys
tl = sys.argv[1]
# Minutes held back for the final checkpoint write and the resubmit below.
MARGIN_MIN = 25.0
# Slurm renders a limit as [days-]HH:MM:SS or, under one hour, MM:SS.
days, _, clock = tl.rpartition("-")
parts = clock.split(":")
try:
    h, m, s = ([0.0] * (3 - len(parts)) + [float(x) for x in parts])
    total = (float(days) if days else 0.0) * 1440.0 + h * 60.0 + m + s / 60.0
except ValueError:
    print(0)                      # UNLIMITED or unparseable: no budget
else:
    print(f"{max(total - MARGIN_MIN, 0.0):.1f}")
PY
)
fi
echo "    wall budget ${wall_stop_min} min (limit ${tl:-?})"
ISMIP7_WALL_STOP_MIN="$wall_stop_min" srun -n "${SLURM_NTASKS}" python -u "$DRIVER"
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
# The chain is built on auto-resume: a successor inherits the environment and
# no ISMIP7_RESTART, so with resume off it would cold-start at t_start, end at
# the same year, and do it again forever. One job, and say why.
if [ "$ISMIP7_AUTO_RESUME" != "1" ]; then
    echo "chain: auto-resume is off (ISMIP7_AUTO_RESUME=$ISMIP7_AUTO_RESUME)," \
         "so a successor would cold-start and repeat these years. Not"
    echo "chain: resubmitting. This run stops here; set ISMIP7_AUTO_RESUME=1" \
         "to chain, or continue it by hand with ISMIP7_RESTART."
    exit 0
fi

# Ask the driver's own final checkpoint, which is both the state the next job
# would resume from and the record of why this one stopped: `t_yr` is the year
# actually reached and `stalled` is the driver's verdict that the solver gave
# up. It must be THIS job's checkpoint: the results directory is flat and
# shared, and the experiment name is built inside the driver, so neither the
# newest file nor a name rebuilt in the shell identifies this run. The driver
# prints "Saved: <final.h5>" on every exit, including an early stop, so take
# the path from our own Slurm log.
read -r LAST_YEAR STALLED START_YEAR <<EOF
$(NM_LOG="logs/ismip7_fwd_${SLURM_JOB_ID}.out" python - <<'PY'
import os, re
log = os.environ["NM_LOG"]
chk_fn, start_year = None, ""
with open(log, errors="replace") as f:
    for line in f:
        if line.startswith("Saved: ") and line.rstrip().endswith("_final.h5"):
            chk_fn = line[len("Saved: "):].strip()
        else:
            m = re.search(r"Restart: evolved geometry .*t_yr=([-\d.eE+]+)", line)
            if m:
                start_year = m.group(1)
last, stalled = "", 0
if chk_fn and os.path.exists(chk_fn):
    import h5py
    with h5py.File(chk_fn, "r") as h:
        t = h["/"].attrs.get("t_yr")
        last = "" if t is None else repr(float(t))
        stalled = int(h["/"].attrs.get("stalled", 0))
print(last or "-", stalled, start_year or "-")
PY
)
EOF
T_END="${ISMIP7_T_END}"
if [ "$STALLED" = "1" ]; then
    echo "chain: the solver stalled (rescue ladder + subcycles exhausted) at ${LAST_YEAR}."
    echo "chain: NOT resubmitting. A stall is not a run that needs more time;"
    echo "chain: resuming would re-attempt the same years and stall again."
    exit 1
fi
if [ "$LAST_YEAR" = "-" ] || [ -z "$LAST_YEAR" ]; then
    echo "chain: could not read this job's final checkpoint from" \
         "logs/ismip7_fwd_${SLURM_JOB_ID}.out; not resubmitting"
    exit 0
fi
# Finished before anything else: a completed run re-submitted under
# ISMIP7_AUTO_RESUME resumes at t_end, runs zero steps and ends at t_end, and
# that is success, not the no-progress case below.
done_yet=$(python -c "print(1 if float('$LAST_YEAR') >= float('$T_END') - 1e-6 else 0)")
if [ "$done_yet" = "1" ]; then
    echo "chain: reached t_end ${T_END} at ${LAST_YEAR}, finished"
    exit 0
fi
# A clean exit short of t_end that also made no progress is not the wall-clock
# case the chain exists for. Allow exactly one retry, in case the job died for
# a transient reason before its first step, then stop rather than loop.
RETRY="${ISMIP7_CHAIN_NOPROGRESS:-0}"
if [ "$START_YEAR" != "-" ]; then
    progressed=$(python -c "print(1 if float('$LAST_YEAR') > float('$START_YEAR') + 1e-6 else 0)")
    if [ "$progressed" = "0" ]; then
        if [ "$RETRY" != "0" ]; then
            echo "chain: no progress again (resumed at ${START_YEAR}, ended at" \
                 "${LAST_YEAR}) after a retry; stopping"
            exit 1
        fi
        echo "chain: no progress (resumed at ${START_YEAR}, ended at ${LAST_YEAR});" \
             "retrying once"
        RETRY=1
    else
        RETRY=0
    fi
fi

echo "chain: at ${LAST_YEAR} of ${T_END}, resubmitting"
# Without the derived budget: the successor may land on a partition with a
# different limit and must derive its own. Without the restart either: an
# explicit one-off ISMIP7_RESTART wins over auto-resume in every driver, so
# inheriting it would send each link back to the same starting checkpoint and
# the chain would never advance past the first job's end year.
unset ISMIP7_WALL_STOP_MIN ISMIP7_RESTART
sbatch --dependency=afterok:"${SLURM_JOB_ID}" \
       --export=ALL,ISMIP7_CHAIN_NOPROGRESS="$RETRY" \
       antarctica/scripts/batch_runners/nots_projection.sbatch
