#!/bin/bash
# Shared environment for ISMIP7 icepack2 batch jobs, on any cluster.
# Sourced by every job script and by submit.sh. Nothing here submits anything.
#
# One file per cluster lives in sites/ and answers three questions: where is
# Firedrake, what does the scheduler want, and where does the data live. This
# file picks the right one, checks it is filled in, and turns it into an
# environment. Rice NOTS, IU Quartz and a UChicago Midway stub ship with the
# repository; sites/template.sh is the blank to copy for anywhere else.
#
#   ISMIP7_SITE=iu_quartz sbatch ...     # name it
#   sbatch ...                           # or let the hostname choose
#
# Every value a site file sets can be overridden per submission, because they
# are all written as ${VAR:-default}: `ISMIP7_PART_LONG=debug submit.sh ...`
# works without editing anything.
set -u

_ISMIP7_BR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ismip7_site_file() {
    local host want f name match pat
    host="$(hostname -f 2>/dev/null || hostname)"
    want="${ISMIP7_SITE:-}"
    if [ -n "$want" ]; then
        f="$_ISMIP7_BR_DIR/sites/$want.sh"
        [ -r "$f" ] && { echo "$f"; return 0; }
        echo "ERROR: ISMIP7_SITE=$want but $f is not readable." >&2
        echo "       Available: $(cd "$_ISMIP7_BR_DIR/sites" && ls *.sh | sed 's/\.sh$//' | tr '\n' ' ')" >&2
        return 2
    fi
    for f in "$_ISMIP7_BR_DIR"/sites/*.sh; do
        name="$(basename "$f" .sh)"
        [ "$name" = "template" ] && continue
        match="$(sed -nE 's/^ISMIP7_SITE_MATCH="([^"]*)".*/\1/p' "$f" | head -1)"
        # The patterns are globs for `case` to match the hostname against, so
        # keep the shell from expanding them against the invoking directory.
        set -f
        for pat in $match; do
            # shellcheck disable=SC2254
            case "$host" in $pat) set +f; echo "$f"; return 0 ;; esac
        done
        set +f
    done
    echo "ERROR: no site definition matches host '$host'." >&2
    echo "       Copy $_ISMIP7_BR_DIR/sites/template.sh to sites/<name>.sh, fill it in," >&2
    echo "       then submit with ISMIP7_SITE=<name> (or add the hostname to its" >&2
    echo "       ISMIP7_SITE_MATCH so it is found automatically)." >&2
    return 2
}

_ISMIP7_SITE_FILE="$(ismip7_site_file)" || exit 2
# shellcheck disable=SC1090
. "$_ISMIP7_SITE_FILE"

# The site file must supply these six. Everything else defaults below.
ISMIP7_REQUIRED="ISMIP7_FIREDRAKE ISMIP7_PART_LONG ISMIP7_PART_SHORT ISMIP7_PART_DEBUG ISMIP7_REPO ISMIP7_WORK"
# build_firedrake_rice.sbatch creates the venv, so it runs where ISMIP7_FIREDRAKE
# has nothing to name yet. It asks for the rest.
ISMIP7_REQUIRED_BUILD="ISMIP7_PART_LONG ISMIP7_PART_SHORT ISMIP7_PART_DEBUG ISMIP7_REPO ISMIP7_WORK"

ismip7_site_require() {
    local missing=""
    local v
    for v in ${1:-$ISMIP7_REQUIRED}; do
        [ -z "${!v:-}" ] && missing="$missing $v"
    done
    if [ -n "$missing" ]; then
        echo "ERROR: site '$ISMIP7_SITE_NAME' ($_ISMIP7_SITE_FILE) is missing:$missing" >&2
        echo "       Fill them in there, or export them for this submission." >&2
        exit 2
    fi
}

# The resource request of the running job, as sbatch flags, for a chain link
# that resubmits its own script. The job scripts carry no resource directives
# (submit.sh composes them from the site file), so a successor submitted from
# inside a job would otherwise land on the cluster's defaults: wrong partition,
# wrong wall limit, and often one task, which the runners refuse outright.
# --hint=nomultithread and the job name are not readable back from scontrol.
# sbatch takes its name default from SBATCH_JOB_NAME, which --export=ALL does
# not carry, so a successor would show up in squeue as the script filename.
# Both are restated here; submit.sh passes them for the first link.
# Both chains call this, so the rule lives here once.
ismip7_chain_resources() {
    local info feat tlim memn
    info="$(scontrol show job "${SLURM_JOB_ID:-}" 2>/dev/null)"
    feat="$(echo "$info" | grep -oE 'Features=[^ ]+' | cut -d= -f2)"
    [ "$feat" = "(null)" ] && feat=""
    tlim="$(echo "$info" | grep -oE 'TimeLimit=[^ ]+' | cut -d= -f2)"
    memn="$(echo "$info" | grep -oE 'MinMemoryNode=[^ ]+' | cut -d= -f2)"
    ISMIP7_CHAIN_RES=(-N "${SLURM_JOB_NUM_NODES:-1}" -n "${SLURM_NTASKS:-1}"
                      --cpus-per-task="${SLURM_CPUS_PER_TASK:-1}"
                      --hint=nomultithread)
    [ -n "${SLURM_JOB_NAME:-}" ] && ISMIP7_CHAIN_RES+=(-J "$SLURM_JOB_NAME")
    [ -n "${SLURM_JOB_PARTITION:-}" ] && ISMIP7_CHAIN_RES+=(-p "$SLURM_JOB_PARTITION")
    [ -n "$feat" ] && ISMIP7_CHAIN_RES+=(-C "$feat")
    [ -n "$tlim" ] && ISMIP7_CHAIN_RES+=(--time="$tlim")
    [ -n "$memn" ] && ISMIP7_CHAIN_RES+=(--mem="$memn")
    [ -n "${SLURM_JOB_ACCOUNT:-}" ] && ISMIP7_CHAIN_RES+=(-A "$SLURM_JOB_ACCOUNT")
    return 0
}

ismip7_load_modules() {
    command -v module >/dev/null 2>&1 || return 0
    module purge 2>/dev/null || true
    if [ -n "${ISMIP7_MODULE_USE:-}" ]; then
        # shellcheck disable=SC2086
        module use ${ISMIP7_MODULE_USE}
    fi
    if [ -n "${ISMIP7_MODULES:-}" ]; then
        # shellcheck disable=SC2086
        module load ${ISMIP7_MODULES}
    fi
}

ismip7_activate() {
    ismip7_site_require
    ismip7_load_modules
    if [ ! -r "$ISMIP7_FIREDRAKE" ]; then
        echo "ERROR: ISMIP7_FIREDRAKE is unreadable: '$ISMIP7_FIREDRAKE'" >&2
        echo "       Set it in $_ISMIP7_SITE_FILE, or export it before submitting." >&2
        exit 2
    fi
    # shellcheck disable=SC1090
    . "$ISMIP7_FIREDRAKE"
    export OMP_NUM_THREADS=1          # one thread per rank; the solver is MPI-parallel
    # Each rank compiles UFL kernels; a shared cache on a networked filesystem
    # corrupts under concurrent writes (seen locally: "undefined symbol:
    # wrap_form0_cell_integral"). Every job gets its own, keyed on the job id,
    # so a chain link killed mid compile cannot hand its successor a truncated
    # object through --export=ALL.
    export PYOP2_CACHE_DIR="${SCRATCH:-$HOME}/.pyop2_cache/${SLURM_JOB_ID:-manual}"
    mkdir -p "$PYOP2_CACHE_DIR"
}

# --- job size ------------------------------------------------------------
# Only the fields ismip7_site_require checks have to come from the site file.
# Everything below carries a working default, so a site file written from that
# list alone still composes a complete submission.
ISMIP7_TASKS="${ISMIP7_TASKS:-16}"
ISMIP7_MEM="${ISMIP7_MEM:-120G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-2-00:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-1-00:00:00}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-}"

# A site that needs one number sets ISMIP7_TASKS/ISMIP7_MEM and both kinds take
# it. A site with measured per-kind values sets the pair. At Rice the forward
# was measured at 12 ranks and the inversion at 32, so running the forward at
# the inversion's size would be an unvalidated rank count on a narrower set of
# nodes.
ISMIP7_TASKS_INV="${ISMIP7_TASKS_INV:-$ISMIP7_TASKS}"
ISMIP7_MEM_INV="${ISMIP7_MEM_INV:-$ISMIP7_MEM}"
ISMIP7_TASKS_FWD="${ISMIP7_TASKS_FWD:-$ISMIP7_TASKS}"
ISMIP7_MEM_FWD="${ISMIP7_MEM_FWD:-$ISMIP7_MEM}"

# The node feature splits the same way. An inversion needs the partition with
# the memory, and the rest of the kinds go wherever the measured timings and
# the -march=native build came from.
ISMIP7_CONSTRAINT_INV="${ISMIP7_CONSTRAINT_INV:-${ISMIP7_CONSTRAINT:-}}"
ISMIP7_CONSTRAINT_FWD="${ISMIP7_CONSTRAINT_FWD:-${ISMIP7_CONSTRAINT:-}}"

# --- repository and data ------------------------------------------------
export ISMIP7_DATA_ROOT="${ISMIP7_DATA_ROOT:-$ISMIP7_REPO/ISMIP7/AIS}"

# --- model configuration shared by every run ----------------------------
# The 2500 m configuration, except for the friction law: every inversion now
# runs regularized Coulomb. Budd's shelf gate was a sign test on the roundoff
# residue of N, so every Budd MAP predating that fix has to be re-inverted;
# set ISMIP7_FRICTION=budd explicitly for those re-inversions.
export ISMIP7_GEOMETRY_SPACE="${ISMIP7_GEOMETRY_SPACE:-dg0}"
export ISMIP7_FRICTION="${ISMIP7_FRICTION:-regularized_coulomb}"
export ISMIP7_N_FLOW="${ISMIP7_N_FLOW:-3.0}"
export ISMIP7_LC="${ISMIP7_LC:-2500}"
export ISMIP7_LC_COARSE="${ISMIP7_LC_COARSE:-64000}"
export ISMIP7_MESH="${ISMIP7_MESH:-$ISMIP7_REPO/antarctica/mesh/antarctica_64000_2500.msh}"
# The MAP this configuration writes and reads: named from the law, so switching
# ISMIP7_FRICTION switches the file and a Budd re-inversion cannot land on the
# RC MAP. One variable for both halves of the workflow - inversion.sbatch
# writes it, projection.sbatch loads it - because the runners' `logvelnet`
# name is not one a forward can derive for itself.
. "$_ISMIP7_BR_DIR/../ismip7_names.sh"
export ISMIP7_MAP_DEFAULT="${ISMIP7_MAP_DEFAULT:-$ISMIP7_REPO/antarctica/mesh/$(ismip7_map_basename "$ISMIP7_FRICTION" "$ISMIP7_LC")}"

ismip7_banner() {
    echo "=== $(date -Is)  job ${SLURM_JOB_ID:-none} on $(hostname) ==="
    echo "    site    $ISMIP7_SITE_NAME ($(basename "$_ISMIP7_SITE_FILE"))"
    echo "    sched   partition ${SLURM_JOB_PARTITION:-?}  nodes ${SLURM_JOB_NUM_NODES:-?}" \
         " ntasks ${SLURM_NTASKS:-?}  mem ${SLURM_MEM_PER_NODE:-?}M" \
         "${ISMIP7_ACCOUNT:+ account $ISMIP7_ACCOUNT}"
    echo "    node    $(scontrol show node "$(hostname -s)" 2>/dev/null \
                      | tr ' ' '\n' | grep -E '^(CPUTot|RealMemory|ActiveFeatures)=' \
                      | paste -sd' ' || echo 'unknown')"
    echo "    repo    $ISMIP7_REPO"
    echo "    mesh    $ISMIP7_MESH"
    echo "    lc=$ISMIP7_LC lc_coarse=$ISMIP7_LC_COARSE geometry=$ISMIP7_GEOMETRY_SPACE" \
         "friction=$ISMIP7_FRICTION n=$ISMIP7_N_FLOW"
}
