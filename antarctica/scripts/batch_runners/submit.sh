#!/bin/bash
# Submit an ISMIP7 job with this cluster's own scheduler settings.
#
#   submit.sh inversion  [options] [KEY=VALUE ...]
#   submit.sh projection [options] [KEY=VALUE ...]
#   submit.sh smoke|probe|verify|build [options] [KEY=VALUE ...]
#
# The resources come from sites/<this cluster>.sh (see site_env.sh), so the
# same command line works at Rice, at IU and anywhere a site file exists. Any
# KEY=VALUE is exported into the job, which is how run settings are chosen:
#
#   submit.sh inversion ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 \
#       ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_5000_2000_buffered0.msh
#   submit.sh projection ISMIP7_EXPERIMENT=ssp585_cesm_waccm ISMIP7_OUTPUT=1
#   ISMIP7_SITE=iu_quartz submit.sh inversion --dry-run
#
# Options (each overrides the site default, in either spelling, --mem 240G or
# --mem=240G):
#   --tasks N --mem 240G --time 1-00:00:00
#   --partition P --constraint C --account A --name JOBNAME
#   --dry-run    print the sbatch command and stop
#
# Why a wrapper rather than #SBATCH lines in each script: a directive is parsed
# before any shell runs, so it cannot read a site file, and a mismatched pair (a
# header asking for 12 tasks per node while the command line asks for 32 in
# total) is rejected with "Requested node configuration is not available". The
# job scripts therefore carry no resource directives at all; this composes them.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
. "$here/site_env.sh"
ismip7_site_require

kind="${1:-}"; shift || true
case "$kind" in
    inversion)  script=inversion.sbatch;       part="$ISMIP7_PART_LONG";  time="$ISMIP7_TIME_INV"; tasks="$ISMIP7_TASKS_INV"; mem="$ISMIP7_MEM_INV"; cons="$ISMIP7_CONSTRAINT_INV"; name=ismip7_inv ;;
    projection) script=projection.sbatch;      part="$ISMIP7_PART_SHORT"; time="$ISMIP7_TIME_FWD"; tasks="$ISMIP7_TASKS_FWD"; mem="$ISMIP7_MEM_FWD"; cons="$ISMIP7_CONSTRAINT_FWD"; name=ismip7_fwd ;;
    smoke)      script=smoke.sbatch;           part="$ISMIP7_PART_DEBUG"; time=00:45:00;           tasks=4;              mem=24G;           cons="$ISMIP7_CONSTRAINT_FWD"; name=ismip7_smoke ;;
    verify)     script=verify.sbatch;          part="$ISMIP7_PART_DEBUG"; time=00:15:00;           tasks=4;              mem=16G;           cons="$ISMIP7_CONSTRAINT_FWD"; name=fd_verify ;;
    probe)      script=partition_probe.sbatch; part="$ISMIP7_PART_DEBUG"; time=01:00:00;           tasks="$ISMIP7_TASKS"; mem=64G;          cons="$ISMIP7_CONSTRAINT_FWD"; name=ismip7_part ;;
    build)      script=build_firedrake.sbatch; part="$ISMIP7_PART_SHORT"; time=12:00:00;           tasks=1;              mem=48G;           cons="$ISMIP7_CONSTRAINT_FWD"; name=fd_build ;;
    *) sed -n '2,${/^#/!q;s/^# \{0,1\}//p;}' "$0"; exit 2 ;;
esac

constraint="$cons"
account="${ISMIP7_ACCOUNT:-}"
dry=0
exports=()
while [ $# -gt 0 ]; do
    # --opt=value is the spelling sbatch itself takes, so split it and let the
    # same branches handle both forms.
    case "$1" in
        --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}" ;;
    esac
    case "$1" in
        --tasks)      tasks="$2"; shift 2 ;;
        --mem)        mem="$2"; shift 2 ;;
        --time)       time="$2"; shift 2 ;;
        --partition)  part="$2"; shift 2 ;;
        --constraint) constraint="$2"; shift 2 ;;
        --account)    account="$2"; shift 2 ;;
        --name)       name="$2"; shift 2 ;;
        --dry-run)    dry=1; shift ;;
        *=*)          exports+=("$1"); shift ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

# The build job asks for whole cores rather than ranks.
cpus_per_task=1
if [ "$kind" = build ]; then cpus_per_task=16; fi

export_list="ALL,ISMIP7_SITE=$ISMIP7_SITE_NAME,ISMIP7_REPO=$ISMIP7_REPO"
for kv in ${exports+"${exports[@]}"}; do export_list="$export_list,$kv"; done

cmd=(sbatch --parsable
     -J "$name"
     -p "$part"
     --nodes=1 --ntasks-per-node="$tasks" --cpus-per-task="$cpus_per_task"
     --hint=nomultithread
     --mem="$mem" --time="$time"
     --export="$export_list"
     "$here/$script")
[ -n "$constraint" ] && cmd=("${cmd[@]:0:1}" -C "$constraint" "${cmd[@]:1}")
[ -n "$account" ] && cmd=("${cmd[@]:0:1}" -A "$account" "${cmd[@]:1}")

printf 'site %s: ' "$ISMIP7_SITE_NAME"; printf '%q ' "${cmd[@]}"; echo
[ "$dry" = 1 ] && exit 0
cd "$ISMIP7_REPO"
mkdir -p logs
"${cmd[@]}"
