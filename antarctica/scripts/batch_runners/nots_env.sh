#!/bin/bash
# Shared environment for ISMIP7 icepack2 jobs on Rice NOTS (deepsC condo).
# Sourced by every nots_*.sbatch script. Nothing here submits anything.
#
# ACCESS AS OF 2026-09-12, all measured on the login node.
#
# Your Slurm association EXISTS and your default account is `commons`, so -A is
# not needed. It grants four partitions. It does NOT grant deepsC.
#
#   partition  wall limit   your cap (per user)      nodes  what it is for
#   commons    1 day        1536 cpu / 2304 GB       142    forwards, 2500 m inversions
#   long       3 days       384 cpu / 1 TB            61    the 2 km inversions
#   debug      30 min       384 cpu / 576 GB, 1 job    2    smoke tests
#   scavenge   1 hour       300 jobs, preemptible    154    probes
#   deepsC     infinite     NO ACCESS                 52    the EEPS condo
#
# deepsC needs the separate `deepsc` Slurm account, held by the EEPS condo
# owner. A --test-only submit there still fails with "Invalid account or
# account/partition combination" while the same submit succeeds on all four
# partitions above. Getting added is a request to the condo owner, not a CRC
# ticket. Until then set NOTS_PARTITION=long or commons.
#
# HARDWARE. deepsC is 52 uniform nodes: 40 physical cores (2x20,
# ThreadsPerCore=2, so Slurm advertises 80), 187135 MB (~182 GiB), features
# `cascadelake,opath`.
#
# AN EARLIER VERSION OF THIS FILE CLAIMED deepsC HAD 768 GB AND 1.5 TB TIERS.
# That was wrong. Every deepsC node is 187 GB. The big-memory Cascade Lake nodes
# live elsewhere, and usefully for us they are in `long`: one 1.5 TB node
# (bc11u27n4, cascadelake), plus 515 GB Sapphire Rapids at 192 cpus. So the
# 2 km / 5 km-interior inversion at ~255 GB runs on `long`, inside the 1 TB
# per-user cap, and could never have run on deepsC.
#
# Cores are hyperthreaded: 40 physical, 80 logical. The solver is memory-
# bandwidth bound, so ranks go on PHYSICAL cores (--hint=nomultithread);
# treating 80 as available would halve per-rank bandwidth.
#
# CONTENTION is real: 624 jobs pending on commons, 216 on long, 182 on deepsC
# when measured. Size jobs to start, not to be perfect.
#
# THREE VALUES MUST BE FILLED IN BEFORE THE FIRST SUBMISSION. They are the
# only site-specific unknowns, and they are deliberately left empty rather
# than guessed, because a wrong module name fails minutes into a queued job.
#
#   NOTS_ACCOUNT     RESOLVED 2026-09-12: the default account `commons` works
#                    and -A is unnecessary. Leave empty unless you are charging
#                    a condo (then NOTS_ACCOUNT=deepsc, once granted).
#   NOTS_FIREDRAKE   the activate script of the Firedrake venv on NOTS, e.g.
#                    /projects/<group>/sw/firedrake-2026.4.1/bin/activate
#   NOTS_MODULES     the module loads that Firedrake venv was BUILT against.
#                    They must match the build exactly or MPI will mismatch.
#
# Everything else below is measured rather than assumed: see
# antarctica/scripts/batch_runners/readme.md for where each number came from.

set -u

NOTS_ACCOUNT="${NOTS_ACCOUNT:-}"          # empty = use the default (commons)

# Partition every job goes to. `long` is the default because it is the only
# one that both fits a 2500 m inversion in its 3-day limit and can host the
# 2 km runs. Override per submission: NOTS_PARTITION=commons for <24 h work,
# debug for a smoke test, deepsC once the condo grants access.
NOTS_PARTITION="${NOTS_PARTITION:-long}"
NOTS_FIREDRAKE="${NOTS_FIREDRAKE:-/projects/ah301/sw/venv-firedrake/bin/activate}"

# Cascade Lake feature name, if deepsC defines one. Empty = do not constrain
# (safe, since every deepsC node is Cascade Lake today). Set it if the condo
# ever gains mixed hardware.
# NOTE: this variable is INFORMATIONAL ONLY. #SBATCH directives are parsed
# before the shell runs, so they cannot read a variable; the constraint is
# written literally as `#SBATCH -C cascadelake` in each job script. Override it
# per submission on the command line, which DOES beat the directive:
#     sbatch -C sapphirerapids ...        # newer, 192 cpus, 257-515 GB
#     sbatch --constraint= ...            # any node, fastest to start
# Kept here only so the banner can report what the scripts pin.
NOTS_CONSTRAINT="${NOTS_CONSTRAINT:-cascadelake}"

# --- modules -----------------------------------------------------------
# Replace this block with the loads used to build the venv. Left empty so a
# misconfigured job fails loudly at the check below rather than silently
# running against the wrong MPI.
# The exact set nots_build_firedrake.sbatch built against (job 1288484,
# 2026-09-12), plus libGLU: the gmsh wheel dlopens libGLU.so.1 at import and
# the compute nodes do not have it. Change the rest only if the venv is rebuilt.
NOTS_MODULES="${NOTS_MODULES:-foss/2023b Python/3.11.5 CMake/3.27.6 M4/1.4.19 flex/2.6.4 Bison/3.8.2 libevent/2.1.12 zstd/1.5.5 Szip/2.1.1 libaec/1.0.6 libGLU/9.0.3}"

nots_load_modules() {
    module purge
    if [ -n "${NOTS_MODULES:-}" ]; then
        # shellcheck disable=SC2086
        module load ${NOTS_MODULES}
    fi
}

nots_activate() {
    nots_load_modules
    if [ -z "$NOTS_FIREDRAKE" ] || [ ! -r "$NOTS_FIREDRAKE" ]; then
        echo "ERROR: NOTS_FIREDRAKE is unset or unreadable: '$NOTS_FIREDRAKE'" >&2
        echo "       Set it in nots_env.sh, or export it before sbatch." >&2
        exit 2
    fi
    # shellcheck disable=SC1090
    . "$NOTS_FIREDRAKE"
    export OMP_NUM_THREADS=1          # one thread per rank; the solver is MPI-parallel
    # Each rank compiles UFL kernels; a shared cache on a networked filesystem
    # corrupts under concurrent writes (seen locally: "undefined symbol:
    # wrap_form0_cell_integral"). Give every job its own.
    export PYOP2_CACHE_DIR="${SCRATCH:-$HOME}/.pyop2_cache/${SLURM_JOB_ID:-manual}"
    mkdir -p "$PYOP2_CACHE_DIR"
}

# --- repository and data ------------------------------------------------
# /projects/ah301 exists and is writable (20 TB share, 63% used). /home is a
# 10 TB NFS export and /scratch is 1.1 PB but 88% full and purged. The forcing
# tree is 316 GB, so it belongs on /projects, not /home.
ISMIP7_REPO="${ISMIP7_REPO:-/projects/ah301/ismip7}"
export ISMIP7_DATA_ROOT="${ISMIP7_DATA_ROOT:-$ISMIP7_REPO/ISMIP7/AIS}"

# --- model configuration shared by every run ----------------------------
# The 2500 m configuration, except for the friction law: every inversion now
# runs regularized Coulomb. Budd's shelf gate was a sign test on the roundoff
# residue of N, so every Budd MAP predating that fix has to be re-inverted;
# set ISMIP7_FRICTION=budd explicitly for those re-inversions.
export ISMIP7_GEOMETRY_SPACE="${ISMIP7_GEOMETRY_SPACE:-dg0}"
export ISMIP7_FRICTION="${ISMIP7_FRICTION:-regularized_coulomb}"
# Filename tag for the law, the same mapping icepack2_tools.naming._FRICTION_TAGS
# uses. Every default MAP name interpolates it, so switching the law switches
# the file the run writes and a Budd re-inversion cannot land on an RC MAP.
case "$ISMIP7_FRICTION" in
    regularized_coulomb) export ISMIP7_FRICTION_TAG=_rc ;;
    budd)                export ISMIP7_FRICTION_TAG=_budd ;;
    *)                   export ISMIP7_FRICTION_TAG= ;;
esac
export ISMIP7_N_FLOW="${ISMIP7_N_FLOW:-3.0}"
export ISMIP7_LC="${ISMIP7_LC:-2500}"
export ISMIP7_LC_COARSE="${ISMIP7_LC_COARSE:-64000}"
export ISMIP7_MESH="${ISMIP7_MESH:-$ISMIP7_REPO/antarctica/mesh/antarctica_64000_2500.msh}"

nots_banner() {
    echo "=== $(date -Is)  job ${SLURM_JOB_ID:-none} on $(hostname) ==="
    echo "    partition ${SLURM_JOB_PARTITION:-?}  nodes ${SLURM_JOB_NUM_NODES:-?}" \
         " ntasks ${SLURM_NTASKS:-?}  mem ${SLURM_MEM_PER_NODE:-?}M"
    echo "    node    $(scontrol show node "$(hostname -s)" 2>/dev/null \
                      | tr ' ' '\n' | grep -E '^(CPUTot|RealMemory|ActiveFeatures)=' \
                      | paste -sd' ' || echo 'unknown')"
    echo "    repo    $ISMIP7_REPO"
    echo "    mesh    $ISMIP7_MESH"
    echo "    lc=$ISMIP7_LC lc_coarse=$ISMIP7_LC_COARSE geometry=$ISMIP7_GEOMETRY_SPACE" \
         "friction=$ISMIP7_FRICTION n=$ISMIP7_N_FLOW"
}
