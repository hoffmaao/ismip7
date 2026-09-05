#!/usr/bin/env bash
# Resource-gated launcher for the 2 km / 5 km-interior Budd inversion.
#
# Mesh: antarctica_5000_2000_buffered0.msh, generated Sep 2026 with the
# current mesh_antarctica.py sizing. 925,183 vertices, 1,835,718 cells,
# element sizes 2000-5000 m. This is the first Antarctic mesh here whose
# INTERIOR is actually resolved: every earlier mesh graded a 2-2.5 km margin
# against a 64 km interior, so interior friction and fluidity were being
# inverted on 64 km cells. It ships its own boundary_ids sidecar (320 calving
# + 319 other markers); the solver hard-errors without a matching one, and no
# previous 2 km mesh had it.
#
# WHY THIS IS GATED. The momentum solve is a direct MUMPS LU of the mixed
# V x Sigma x T system, and there is no iterative fallback anywhere in the
# inversion. Measured on this machine at 2500 m (Z = 3,470,017 dof, 12 ranks):
# 28.4 GB summed over ranks, 1181 s of setup, ~840 s per L-BFGS iteration.
# This mesh is Z = 11,028,956 dof, 3.18x that, so expect roughly 90-110 GB and
# 4-5 days for a converged run. On 2026-09-04 systemd-oomd killed a whole
# terminal (227 processes, a 23-hour inversion with it) when the user slice hit
# 91.58% memory pressure; the gia-icepack fleet alone was holding 255 GB. So
# this must NOT be launched into a busy machine. The gate below waits for real
# headroom instead.
#
#   setsid nohup antarctica/scripts/launch_inversion_2km_when_ready.sh >/dev/null 2>&1 &
# DRYRUN=1 prints one resource reading and exits.
# rm the lock under antarctica/results/logs/ to re-arm.
set -u

REPO="${REPO:-/media/andrew/wd1/projects/ismip7}"
PY="${PY:-$HOME/venv-firedrake-2026/bin/python}"
NRANKS="${NRANKS:-12}"          # matches the 2500 m memory measurement
NEED_GB="${NEED_GB:-160}"       # ~110 GB for the factorisation + headroom
NEED_CORES="${NEED_CORES:-20}"
CHECKS="${CHECKS:-5}"
INTERVAL="${INTERVAL:-180}"

LOGDIR="$REPO/antarctica/results/logs"
LOCK="$LOGDIR/inv2km.lock"
LOG="$LOGDIR/inv2km_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$LOGDIR"

if [ -e "$LOCK" ]; then echo "already armed/ran (lock: $LOCK)"; exit 0; fi

read_free_gb() { awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo; }
read_idle_cores() {
  local ncpu load
  ncpu=$(nproc)
  load=$(awk '{print int($1)}' /proc/loadavg)
  echo $(( ncpu - load < 0 ? 0 : ncpu - load ))
}

if [ "${DRYRUN:-0}" = "1" ]; then
  echo "free=$(read_free_gb) GB (need $NEED_GB), idle cores=$(read_idle_cores) (need $NEED_CORES)"
  exit 0
fi

echo "$$" > "$LOCK"
ok=0
while :; do
  if [ "$(read_free_gb)" -ge "$NEED_GB" ] && [ "$(read_idle_cores)" -ge "$NEED_CORES" ]; then
    ok=$((ok + 1))
  else
    ok=0
  fi
  [ "$ok" -ge "$CHECKS" ] && break
  sleep "$INTERVAL"
done

cd "$REPO"
# Settings mirror the converged 2500 m log-velocity inversion
# (inversion_icepack2_budd_n3_dg0_logvel_2500.h5): same misfit, same priors,
# same friction. Only the mesh changes, so the two are comparable and
# score_map.py can put them side by side.
OMP_NUM_THREADS=1 \
ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_BUFFER_M=0 \
ISMIP7_MESH="$REPO/antarctica/mesh/antarctica_5000_2000_buffered0.msh" \
ISMIP7_FRICTION=budd ISMIP7_N_FLOW=3.0 ISMIP7_GEOMETRY_SPACE=dg0 \
ISMIP7_MISFIT_NORM=sigma ISMIP7_SIGMA_U_FLOOR=3 \
ISMIP7_LOG_VEL_WEIGHT=auto \
ISMIP7_DHDT_WEIGHT=1 ISMIP7_GAMMA_THETA=1e5 ISMIP7_GAMMA_PHI=1e5 \
ISMIP7_MAXITER="${ISMIP7_MAXITER:-200}" \
ISMIP7_MAP_OUT="$REPO/antarctica/mesh/inversion_icepack2_budd_n3_dg0_logvel_2000.h5" \
nice -n 5 mpiexec --mca btl self,vader --mca btl_base_warn_component_unused 0 \
  -n "$NRANKS" "$PY" -u antarctica/scripts/inversion_icepack2.py >> "$LOG" 2>&1
echo "INVERSION-2KM-DONE rc=$? $(date)" >> "$LOG"
