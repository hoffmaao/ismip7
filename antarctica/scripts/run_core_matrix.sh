#!/bin/bash
# run_core_matrix.sh - run the ISMIP7 AIS core-experiment matrix (cores 1-11)
# end to end, in protocol dependency order, with load gating and wall retry.
#
# Protocol order (cores 3-10 all branch from the historical endpoint, so the
# historicals must finish first; see experiment.py and control/run.py):
#
#   1,2   historical   CESM2-WACCM / MRI-ESM2-0        1850-2014
#   9,10  CTRL2015     both ESMs, constant 2015 climate   ->2300
#   3,4   ssp370       both ESMs                          ->2100
#   5,6   ssp126       both ESMs                          ->2300
#   7,8   ssp585       both ESMs                          ->2300
#   11    OCX          obs-constrained, cold start     1990-2025
#
# WALL RETRY: the split-step diagnostic solve can exhaust its in-run rescue
# ladder late in a long projection and save-and-stop short of the target year.
# Relaunching from the saved state clears it - a fresh process re-runs the
# n=1->n continuation at the loaded geometry, which the in-run ladder cannot do
# (observed 3/3: ssp585-CESM 2096.7, CTRL-CESM 2268, CTRL-MRI 2250 all resumed
# and two reached 2300). This script therefore relaunches a short run from its
# own final.h5, and keeps doing so while each attempt makes progress, up to
# MAX_ATTEMPTS. A run that stops advancing is left alone and reported.
#
# Usage:
#   antarctica/scripts/run_core_matrix.sh                # full matrix
#   CORES=1,2,9 antarctica/scripts/run_core_matrix.sh    # a subset
#   ISMIP7_RUN_TAG=n3 NRANKS=8 antarctica/scripts/run_core_matrix.sh
#
# Env: ISMIP7_LC (32000), ISMIP7_N_FLOW (3), ISMIP7_RUN_TAG, ISMIP7_DT (0.1),
#      ISMIP7_FRICTION (budd), NRANKS (8), MAX_LOAD (cores-8), MAX_ATTEMPTS (6),
#      CORES (comma list, default all), VENV (~/venv-firedrake-2026).
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO"
R=antarctica/results
mkdir -p "$R/logs"

VENV="${VENV:-$HOME/venv-firedrake-2026}"
NRANKS="${NRANKS:-8}"          # 8 is MUMPS-safe at 32 km; raise for finer meshes
MAX_ATTEMPTS="${MAX_ATTEMPTS:-6}"
NCPU=$(python3 -c 'import os;print(os.cpu_count())')
MAX_LOAD="${MAX_LOAD:-$((NCPU - 8))}"
CORES="${CORES:-1,2,9,10,3,4,5,6,7,8,11}"

export OMP_NUM_THREADS=1
export ISMIP7_LC="${ISMIP7_LC:-32000}"
export ISMIP7_N_FLOW="${ISMIP7_N_FLOW:-3}"
export ISMIP7_FRICTION="${ISMIP7_FRICTION:-budd}"
export ISMIP7_DT="${ISMIP7_DT:-0.1}"
export ISMIP7_OUTPUT_INTERVAL="${ISMIP7_OUTPUT_INTERVAL:-10}"
export ISMIP7_APPARENT_MB="${ISMIP7_APPARENT_MB:-1}"
TAG="${ISMIP7_RUN_TAG:-}"
[ -n "$TAG" ] && export ISMIP7_RUN_TAG="$TAG"
SFX=${TAG:+_$TAG}
LC=$ISMIP7_LC
TS=$(date +%Y%m%d_%H%M%S)
PY="$VENV/bin/python"

# core | label | driver | ESM | experiment-name stem | target year
CORE_SPEC=(
  "1|hist_cesm|historical/cesm_waccm.py|CESM2-WACCM|hist_cesm2_waccm|2014"
  "2|hist_mri|historical/mri_esm2.py|MRI-ESM2-0|hist_mri_esm2_0|2014"
  "9|ctrl_cesm|control/run.py|CESM2-WACCM|ctrl2015_cesm2_waccm|2300"
  "10|ctrl_mri|control/run.py|MRI-ESM2-0|ctrl2015_mri_esm2_0|2300"
  "3|ssp370_cesm|projections/ssp370_cesm_waccm.py|CESM2-WACCM|ssp370_cesm2_waccm|2100"
  "4|ssp370_mri|projections/ssp370_mri_esm2.py|MRI-ESM2-0|ssp370_mri_esm2_0|2100"
  "5|ssp126_cesm|projections/ssp126_cesm_waccm.py|CESM2-WACCM|ssp126_cesm2_waccm|2300"
  "6|ssp126_mri|projections/ssp126_mri_esm2.py|MRI-ESM2-0|ssp126_mri_esm2_0|2300"
  "7|ssp585_cesm|projections/ssp585_cesm_waccm.py|CESM2-WACCM|ssp585_cesm2_waccm|2300"
  "8|ssp585_mri|projections/ssp585_mri_esm2.py|MRI-ESM2-0|ssp585_mri_esm2_0|2300"
  "11|ocx|projections/ocx.py|CESM2-WACCM|ocx|2025"
)

wanted () { [[ ",$CORES," == *",$1,"* ]]; }
csv_of () { echo "$R/${1}${SFX}_${LC}_timeseries.csv"; }
h5_of ()  { echo "$R/${1}${SFX}_${LC}_final.h5"; }

last_year () {  # last simulated year of a timeseries, or empty
  local c="$1"
  [ -f "$c" ] || return 1
  awk -F, 'NR>1 && $1+0>0 {y=$1} END{if(y!="") print y}' "$c"
}

wait_for_load () {
  for _ in $(seq 1 720); do
    local l; l=$(python3 -c 'import os;print(int(os.getloadavg()[0]))')
    [ "$l" -le "$MAX_LOAD" ] && return 0
    sleep 30
  done
  echo "    (load never dropped below $MAX_LOAD; proceeding anyway)"
}

run_core () {  # core label driver esm stem target
  local core="$1" label="$2" driver="$3" esm="$4" stem="$5" target="$6"
  local csv h5 prev now attempt log
  csv=$(csv_of "$stem"); h5=$(h5_of "$stem")

  now=$(last_year "$csv" || echo "")
  if [ -n "$now" ] && awk "BEGIN{exit !($now >= $target)}"; then
    echo "[core $core] $label already at $now (target $target) - skipping"
    return 0
  fi

  for attempt in $(seq 1 "$MAX_ATTEMPTS"); do
    prev="$now"
    wait_for_load
    log="$R/logs/matrix_${TS}_core${core}_${label}_a${attempt}.log"
    # Attempt 1 uses the driver's own restart logic (projections and the CTRL
    # branch from the historical endpoint). Later attempts are wall retries and
    # must resume from this run's own saved state.
    local restart_env=()
    if [ "$attempt" -gt 1 ]; then
      [ -f "$h5" ] || { echo "[core $core] no $h5 to retry from - stopping"; return 1; }
      restart_env=(ISMIP7_RESTART="$REPO/$h5")
      echo "[core $core] $label wall retry $((attempt-1)) from $now"
    else
      echo "[core $core] $label starting (target $target)"
    fi
    env "${restart_env[@]}" ISMIP7_ESM="$esm" \
      mpiexec -n "$NRANKS" "$PY" "antarctica/scripts/$driver" \
      ${TAG:+--tag "$TAG"} > "$log" 2>&1
    now=$(last_year "$csv" || echo "")
    echo "    -> reached ${now:-nothing}  (log ${log##*/})"
    [ -z "$now" ] && { echo "[core $core] produced no output - see log"; return 1; }
    awk "BEGIN{exit !($now >= $target)}" && { echo "[core $core] $label COMPLETE at $now"; return 0; }
    # only retry while the wall is actually moving forward
    if [ -n "$prev" ] && awk "BEGIN{exit !($now <= $prev + 0.001)}"; then
      echo "[core $core] $label stalled at $now (no progress) - giving up"
      return 1
    fi
  done
  echo "[core $core] $label short of target after $MAX_ATTEMPTS attempts (at $now)"
  return 1
}

echo "ISMIP7 core matrix | lc=$LC n=$ISMIP7_N_FLOW tag='${TAG:-none}' dt=$ISMIP7_DT ranks=$NRANKS"
echo "cores: $CORES | max load $MAX_LOAD/$NCPU | wall retries up to $MAX_ATTEMPTS"
echo

for spec in "${CORE_SPEC[@]}"; do
  IFS='|' read -r core label driver esm stem target <<< "$spec"
  wanted "$core" || continue
  run_core "$core" "$label" "$driver" "$esm" "$stem" "$target"
done

echo
echo "=== matrix summary ==="
for spec in "${CORE_SPEC[@]}"; do
  IFS='|' read -r core label driver esm stem target <<< "$spec"
  wanted "$core" || continue
  c=$(csv_of "$stem"); y=$(last_year "$c" || echo "-")
  printf "  core %-2s %-14s %8s / %s\n" "$core" "$label" "$y" "$target"
done

echo
echo "=== observational audit ==="
for spec in "${CORE_SPEC[@]}"; do
  IFS='|' read -r core label driver esm stem target <<< "$spec"
  wanted "$core" || continue
  c=$(csv_of "$stem"); [ -f "$c" ] || continue
  printf "  %-14s %s\n" "$label" "$("$PY" antarctica/scripts/check_ismip6_track.py "$c" 2>&1 | tail -1)"
done

echo
echo "=== ISMIP6 ensemble comparison (projection - same-ESM CTRL) ==="
for spec in "5|ssp126_cesm2_waccm|ctrl2015_cesm2_waccm" "6|ssp126_mri_esm2_0|ctrl2015_mri_esm2_0" \
            "3|ssp370_cesm2_waccm|ctrl2015_cesm2_waccm" "4|ssp370_mri_esm2_0|ctrl2015_mri_esm2_0" \
            "7|ssp585_cesm2_waccm|ctrl2015_cesm2_waccm" "8|ssp585_mri_esm2_0|ctrl2015_mri_esm2_0"; do
  IFS='|' read -r core p c <<< "$spec"
  wanted "$core" || continue
  pc=$(csv_of "$p"); cc=$(csv_of "$c")
  [ -f "$pc" ] && [ -f "$cc" ] || continue
  echo "--- core $core: $p ---"
  "$PY" antarctica/scripts/compare_ismip6.py "$pc" "$cc" 2>&1 | tail -3
done
