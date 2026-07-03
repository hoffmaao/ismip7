#!/usr/bin/env bash
# Resource-gated launcher for the 500 m regularized-Coulomb CTRL2015 run.
#
# Runs antarctica/scripts/control/run.py with ISMIP7_FRICTION=regularized_coulomb
# against inversion_icepack2_rc_500.h5: constant 2015-2029 SMB climatology,
# OI-climatology ocean + per-basin calibrated K, dt=0.1, 2015->T_END.
# This is the clamp-free verification of the RC + h_clamp=0 pipeline: the
# initial state is the true BedMachine geometry and no thickness floor is
# ever applied, so total mass should track (SMB - melt - discharge) with no
# buffer-cell resurrection artifact.
#
# Polls free RAM and idle cores; once BOTH clear their thresholds for STABLE
# consecutive checks, launches detached (nohup mpiexec) and exits. Single-
# instance via a lockfile retained after launch ("already fired" marker —
# remove it to re-arm). Logs live under antarctica/results/logs/ so they
# survive reboots (the Jun 20 RC inversion log died with /tmp).
#
# Run detached so it survives the shell/session:
#   setsid nohup antarctica/scripts/launch_rc_ctrl_when_ready.sh >/dev/null 2>&1 &
# DRYRUN=1 prints one resource reading and exits.
set -u

# ---- configuration (all env-overridable) ----
REPO="${REPO:-/home/andrew/projects/ismip7}"
PY="${PY:-$HOME/venv-firedrake/bin/python}"
# The 500 m aniso mesh does not match the antarctica_<lc_coarse>_<lc>.msh
# default pattern, so pass it explicitly (must be the RC inversion's mesh;
# simulation.py verifies node ordering against the checkpoint).
MESH="${ISMIP7_MESH:-$REPO/antarctica/mesh/antarctica_64000_500_aniso.msh}"
LC="${ISMIP7_LC:-500}"
T_END="${ISMIP7_T_END:-2025}"              # 10-yr verification by default
DT="${ISMIP7_DT:-0.1}"
K_NPZ="${ISMIP7_K_PER_BASIN_NPZ:-$REPO/antarctica/results/calibrated_K_per_basin_2500.npz}"
NRANKS="${NRANKS:-24}"
MIN_FREE_GB="${MIN_FREE_GB:-128}"          # ~7M-DOF MUMPS, conservative
MIN_FREE_CORES="${MIN_FREE_CORES:-28}"     # >= NRANKS with headroom
STABLE="${STABLE:-5}"                      # consecutive passing checks
INTERVAL="${INTERVAL:-180}"                # seconds between checks

LOGDIR="$REPO/antarctica/results/logs"
mkdir -p "$LOGDIR"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="${LOG:-$LOGDIR/rc_ctrl500_${STAMP}.log}"
GATELOG="${GATELOG:-$LOGDIR/rc_ctrl500_gate.log}"
LOCK="${LOCK:-$LOGDIR/rc_ctrl500.lock}"

log(){ echo "[$(date '+%F %T')] $*" | tee -a "$GATELOG"; }

read_avail_gb(){ awk '/MemAvailable/{printf "%d", $2/1024/1024}' /proc/meminfo; }
read_idle_cores(){ awk -v n="$(nproc)" '{c=n-int($1+0.999); print (c<0)?0:c}' /proc/loadavg; }

if [ "${DRYRUN:-}" = "1" ]; then
  echo "mesh=$(basename "$MESH")  ranks=$NRANKS  need: RAM>=${MIN_FREE_GB}G cores>=${MIN_FREE_CORES}"
  echo "now: MemAvailable=$(read_avail_gb)G  idle_cores=$(read_idle_cores)  load1=$(awk '{print $1}' /proc/loadavg)"
  exit 0
fi

# single-instance guard
if [ -e "$LOCK" ]; then
  log "lockfile $LOCK present (gate already running or run already launched). Exiting. rm it to re-arm."
  exit 0
fi
echo "$$" > "$LOCK"
trap 'rm -f "$LOCK"' EXIT

[ -f "$MESH" ]  || { log "ERROR: mesh not found: $MESH"; exit 1; }
[ -x "$PY" ]    || { log "ERROR: python not found: $PY"; exit 1; }
[ -f "$K_NPZ" ] || { log "ERROR: per-basin K npz not found: $K_NPZ"; exit 1; }
[ -f "$REPO/antarctica/mesh/inversion_icepack2_rc_${LC}.h5" ] \
  || { log "ERROR: RC MAP not found: inversion_icepack2_rc_${LC}.h5"; exit 1; }
log "ARMED: RC CTRL2015 lc=$LC dt=$DT t_end=$T_END ranks=$NRANKS mesh=$(basename "$MESH")"
log "  gate: RAM>=${MIN_FREE_GB}G AND idle_cores>=${MIN_FREE_CORES} for ${STABLE} checks @ ${INTERVAL}s; run log -> $LOG"

ok=0
while :; do
  gb="$(read_avail_gb)"; idle="$(read_idle_cores)"; load1="$(awk '{print $1}' /proc/loadavg)"
  if [ "$gb" -ge "$MIN_FREE_GB" ] && [ "$idle" -ge "$MIN_FREE_CORES" ]; then
    ok=$((ok+1))
    log "pass ${ok}/${STABLE}  (avail ${gb}G, idle ${idle} cores, load1 ${load1})"
    if [ "$ok" -ge "$STABLE" ]; then
      log "LAUNCHING RC CTRL -> $LOG"
      cd "$REPO" || { log "ERROR: cd $REPO failed"; exit 1; }
      export OMP_NUM_THREADS=1 ISMIP7_LC="$LC" ISMIP7_MESH="$MESH" \
             ISMIP7_FRICTION=regularized_coulomb \
             ISMIP7_DT="$DT" ISMIP7_T_END="$T_END" \
             ISMIP7_K_PER_BASIN_NPZ="$K_NPZ"
      nohup mpiexec -n "$NRANKS" "$PY" antarctica/scripts/control/run.py > "$LOG" 2>&1 &
      pid=$!
      echo "$pid" >> "$LOCK"
      log "launched: pid=$pid  (lock retained as run marker; rm $LOCK to re-arm)"
      trap - EXIT          # keep the lock so a re-run won't double-launch
      exit 0
    fi
  else
    [ "$ok" -gt 0 ] && log "reset (avail ${gb}G, idle ${idle} cores, load1 ${load1})"
    ok=0
  fi
  sleep "$INTERVAL"
done
