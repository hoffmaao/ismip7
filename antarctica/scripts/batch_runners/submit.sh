#!/bin/bash
# Submit an ISMIP7 job with this cluster's own scheduler settings.
#
#   submit.sh inversion  [options] [KEY=VALUE ...]
#   submit.sh projection [options] [KEY=VALUE ...]
#   submit.sh smoke|probe|verify|build [options] [KEY=VALUE ...]
#   submit.sh script PATH [options] [KEY=VALUE ...]
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
# `script` submits any job script that sources site_core.sh, with this site's
# account, node feature and extra flags. It is how antarctica/Makefile and
# manage_timing_campaign.py submit, so the timing campaign is the same command
# everywhere too. It takes these as well:
#   --queue short|long|debug   the site's partition of that class (short)
#   --cd DIR        submit from ISMIP7_REPO/DIR; PATH is then relative to it
#   --dependency D  passed to sbatch as --dependency=D
#   --wait          block until the job ends and exit with its status
# Standard output is sbatch's own (the job id under --parsable); the composed
# command goes to standard error. A request that one node of this site cannot
# hold (ISMIP7_CORES_PER_NODE, ISMIP7_MEM_PER_NODE) is refused with exit status
# 3 and a `reason=` on standard error. Naming a --constraint yourself skips
# that check, since the limits describe ISMIP7_CONSTRAINT_TIMING's nodes.
#
# Why a wrapper rather than #SBATCH lines in each script: a directive is parsed
# before any shell runs, so it cannot read a site file, and a mismatched pair (a
# header asking for 12 tasks per node while the command line asks for 32 in
# total) is rejected with "Requested node configuration is not available". The
# job scripts therefore carry no resource directives at all; this composes them.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
# Read the site file in a subshell and import only the scheduler fields. The
# model defaults site_env.sh exports (ISMIP7_LC, ISMIP7_MESH, ISMIP7_FRICTION,
# ISMIP7_MAP_DEFAULT and the rest) must stay out of this shell, because
# --export=ALL below would carry them into the job and kill every ${VAR:-...}
# fallback the job scripts resolve for themselves. A KEY=VALUE argument still
# reaches the job, and so does anything the operator exported by hand.
# The build job creates the venv, so it is not asked to name one. It is also
# the one Rice-specific script left, kept as a worked example.
# Nor is a dry run, which starts nothing: it has to work on a laptop.
require_list=ISMIP7_REQUIRED
[ "${1:-}" = build ] && require_list=ISMIP7_REQUIRED_BUILD
for arg in "$@"; do
    [ "$arg" = --dry-run ] && require_list=ISMIP7_REQUIRED_BUILD
done
site_fields=$(
    # shellcheck disable=SC1091
    . "$here/site_env.sh"
    ismip7_site_require "${!require_list}"
    declare -p ISMIP7_SITE_NAME ISMIP7_REPO ISMIP7_ACCOUNT \
               ISMIP7_PART_LONG ISMIP7_PART_SHORT ISMIP7_PART_DEBUG \
               ISMIP7_TASKS ISMIP7_TASKS_INV ISMIP7_TASKS_FWD \
               ISMIP7_MEM_INV ISMIP7_MEM_FWD \
               ISMIP7_TIME_INV ISMIP7_TIME_FWD \
               ISMIP7_CONSTRAINT_INV ISMIP7_CONSTRAINT_FWD \
               ISMIP7_CONSTRAINT_TIMING ISMIP7_SBATCH_EXTRA \
               ISMIP7_CORES_PER_NODE ISMIP7_MEM_PER_NODE
) || exit 2
eval "$site_fields"

kind="${1:-}"; shift || true
case "$kind" in
    inversion)  script=inversion.sbatch;       part="$ISMIP7_PART_LONG";  time="$ISMIP7_TIME_INV"; tasks="$ISMIP7_TASKS_INV"; mem="$ISMIP7_MEM_INV"; cons="$ISMIP7_CONSTRAINT_INV"; name=ismip7_inv ;;
    projection) script=projection.sbatch;      part="$ISMIP7_PART_SHORT"; time="$ISMIP7_TIME_FWD"; tasks="$ISMIP7_TASKS_FWD"; mem="$ISMIP7_MEM_FWD"; cons="$ISMIP7_CONSTRAINT_FWD"; name=ismip7_fwd ;;
    smoke)      script=smoke.sbatch;           part="$ISMIP7_PART_DEBUG"; time=00:45:00;           tasks=4;              mem=24G;           cons="$ISMIP7_CONSTRAINT_FWD"; name=ismip7_smoke ;;
    verify)     script=verify.sbatch;          part="$ISMIP7_PART_DEBUG"; time=00:15:00;           tasks=4;              mem=16G;           cons="$ISMIP7_CONSTRAINT_FWD"; name=fd_verify ;;
    probe)      script=partition_probe.sbatch; part="$ISMIP7_PART_DEBUG"; time=01:00:00;           tasks="$ISMIP7_TASKS"; mem=64G;          cons="$ISMIP7_CONSTRAINT_FWD"; name=ismip7_part ;;
    build)      script=build_firedrake_rice.sbatch; part="$ISMIP7_PART_SHORT"; time=12:00:00;           tasks=1;              mem=48G;           cons="$ISMIP7_CONSTRAINT_FWD"; name=fd_build ;;
    script)
        script="${1:-}"; shift || true
        case "$script" in
            ""|--*|*=*) echo "submit.sh script: name the job script first" >&2; exit 2 ;;
        esac
        part="$ISMIP7_PART_SHORT"; time="$ISMIP7_TIME_FWD"; tasks="$ISMIP7_TASKS_FWD"; mem="$ISMIP7_MEM_FWD"
        cons="$ISMIP7_CONSTRAINT_TIMING"
        name="$(basename "$script")"; name="${name%.*}" ;;
    *) sed -n '2,${/^#/!q;s/^# \{0,1\}//p;}' "$0"; exit 2 ;;
esac

constraint="$cons"
account="$ISMIP7_ACCOUNT"
dry=0
exports=()
subdir=""
dependency=""
wait_flag=0
capped=1
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
        --constraint) constraint="$2"; capped=0; shift 2 ;;
        --account)    account="$2"; shift 2 ;;
        --name)       name="$2"; shift 2 ;;
        --dry-run)    dry=1; shift ;;
        --queue|--cd|--dependency|--wait)
            if [ "$kind" != script ]; then
                echo "$1 belongs to \`submit.sh script\`" >&2; exit 2
            fi
            case "$1" in
                --queue)
                    case "$2" in
                        short) part="$ISMIP7_PART_SHORT" ;;
                        long)  part="$ISMIP7_PART_LONG" ;;
                        debug) part="$ISMIP7_PART_DEBUG" ;;
                        *) echo "--queue takes short, long or debug, not '$2'" >&2; exit 2 ;;
                    esac
                    shift 2 ;;
                --cd)         subdir="$2"; shift 2 ;;
                --dependency) dependency="$2"; shift 2 ;;
                --wait)       wait_flag=1; shift ;;
            esac ;;
        *=*)
            # The key has to be a variable name for the job to read it, and
            # the `env` below would take a leading dash for its own option.
            case "${1%%=*}" in
                ""|[0-9]*|*[![:alnum:]_]*)
                    echo "not a variable name: '${1%%=*}' in $1" >&2; exit 2 ;;
            esac
            exports+=("$1"); shift ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

# The build job asks for whole cores rather than ranks.
cpus_per_task=1
if [ "$kind" = build ]; then cpus_per_task=16; fi

# Submit the script out of the checkout it will run in. The job cds to
# SLURM_SUBMIT_DIR and every path it then resolves, site_env.sh, the driver
# and the chain resubmits, is relative to ISMIP7_REPO, so submitting this
# checkout's copy instead would run link 1 from one tree and links 2..N from
# another. A wrong ISMIP7_REPO now stops here.
#
# `script` names its own path, and with --cd really changes directory rather
# than asking sbatch to (--chdir): the timing scripts cd to SLURM_SUBMIT_DIR
# and resolve scripts/, mesh/ and their %j logs from antarctica/, and --chdir
# leaves SLURM_SUBMIT_DIR where the submission was typed.
script_rel="antarctica/scripts/batch_runners/$script"
[ "$kind" = script ] && script_rel="$script"
if ! cd "$ISMIP7_REPO/$subdir" 2>/dev/null || [ ! -f "$script_rel" ]; then
    echo "ERROR: $ISMIP7_REPO/${subdir:+$subdir/}$script_rel does not exist." >&2
    echo "       ISMIP7_REPO must name the checkout to run; set it in the site" >&2
    echo "       file or export it for this submission." >&2
    exit 2
fi

# sbatch splits --export on commas and documents no quoting, so a list entry
# ISMIP7_SUBCYCLES=1,4,16,64 would reach the job as ISMIP7_SUBCYCLES=1. A value
# holding a comma is set in sbatch's own environment instead, over any value
# the calling shell exported, and ALL carries it whole into the job and its
# chain resubmits.
export_list="ALL,ISMIP7_SITE=$ISMIP7_SITE_NAME,ISMIP7_REPO=$ISMIP7_REPO"
comma_values=()
for kv in ${exports+"${exports[@]}"}; do
    case "$kv" in
        *,*) comma_values+=("$kv") ;;
        *)   export_list="$export_list,$kv" ;;
    esac
done

cmd=(sbatch --parsable
     -J "$name"
     -p "$part"
     --nodes=1 --ntasks-per-node="$tasks" --cpus-per-task="$cpus_per_task"
     --hint=nomultithread
     --mem="$mem" --time="$time"
     --export="$export_list"
     "$script_rel")
[ -n "$constraint" ] && cmd=("${cmd[@]:0:1}" -C "$constraint" "${cmd[@]:1}")
[ -n "$account" ] && cmd=("${cmd[@]:0:1}" -A "$account" "${cmd[@]:1}")
# Flags that go in front of the script path, which has to stay last.
tail_flags=()
# shellcheck disable=SC2206
[ -n "$ISMIP7_SBATCH_EXTRA" ] && tail_flags+=($ISMIP7_SBATCH_EXTRA)
[ -n "$dependency" ] && tail_flags+=(--dependency="$dependency")
[ "$wait_flag" = 1 ] && tail_flags+=(--wait)
last=$((${#cmd[@]} - 1))
cmd=("${cmd[@]:0:$last}" ${tail_flags+"${tail_flags[@]}"} "${cmd[$last]}")
# Put in front last, since the insertions above count positions from sbatch.
[ -n "${comma_values+x}" ] && cmd=(env "${comma_values[@]}" "${cmd[@]}")

# Memory as sbatch spells it (240G, 187000M, a bare number of megabytes), in
# megabytes. Anything else prints nothing and the comparison is skipped.
mem_mb() {
    case "$1" in
        *[!0-9KMGTkmgt]*|"") return 0 ;;
    esac
    local n="${1%[KMGTkmgt]}"
    case "$n" in ""|*[!0-9]*) return 0 ;; esac
    case "$1" in
        *[Kk]) echo $((n / 1024)) ;;
        *[Gg]) echo $((n * 1024)) ;;
        *[Tt]) echo $((n * 1024 * 1024)) ;;
        *)     echo "$n" ;;
    esac
}

if [ "$kind" = script ] && [ "$capped" = 1 ]; then
    reason=""
    if [ -n "$ISMIP7_CORES_PER_NODE" ] && [ "$tasks" -gt "$ISMIP7_CORES_PER_NODE" ]; then
        reason="exceeds_site_cores requested=$tasks limit=$ISMIP7_CORES_PER_NODE"
    fi
    want_mb="$(mem_mb "$mem")"; have_mb="$(mem_mb "$ISMIP7_MEM_PER_NODE")"
    if [ -z "$reason" ] && [ -n "$want_mb" ] && [ -n "$have_mb" ] && [ "$want_mb" -gt "$have_mb" ]; then
        reason="exceeds_site_mem requested=$mem limit=$ISMIP7_MEM_PER_NODE"
    fi
    if [ -n "$reason" ]; then
        echo "site $ISMIP7_SITE_NAME: NOT RUNNABLE $name reason=$reason" >&2
        exit 3
    fi
fi

if [ "$kind" = build ] && [ "$ISMIP7_SITE_NAME" != rice_nots ]; then
    echo "WARNING: $script is the Rice recipe. Its module names are Rice's," >&2
    echo "         so at site '$ISMIP7_SITE_NAME' the job will fail on" >&2
    echo "         'module load'. Copy it to build_firedrake_$ISMIP7_SITE_NAME.sbatch," >&2
    echo "         substitute this cluster's module stack, and submit that." >&2
fi
# A caller of `script` reads the job id from standard output, so the composed
# command goes to standard error there.
if [ "$kind" = script ]; then
    { printf 'site %s: ' "$ISMIP7_SITE_NAME"; printf '%q ' "${cmd[@]}"; echo; } >&2
else
    printf 'site %s: ' "$ISMIP7_SITE_NAME"; printf '%q ' "${cmd[@]}"; echo
fi
[ "$dry" = 1 ] && exit 0
# The runners log under logs/; a `script` job names its own log files.
[ "$kind" = script ] || mkdir -p logs
# Last, so that under --wait this exits with the job's own status.
"${cmd[@]}"
