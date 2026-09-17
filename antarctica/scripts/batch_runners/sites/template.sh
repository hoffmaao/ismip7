# Site definition template. Copy to sites/<your-site>.sh, fill in, and either
# name it in ISMIP7_SITE or let site_core.sh match it by hostname.
#
# A site file describes a CLUSTER and is tracked. What is yours alone on that
# cluster (account, a private build, a mail address) goes in sites/local.env,
# which git ignores; sites/local.env.example is its blank.
#
# Nothing here submits anything: these are the answers to "where is Firedrake,
# what does the scheduler want, and where does the data live" for one cluster.
# site_env.sh refuses to run a job while a REQUIRED value is empty, because a
# wrong module name or a missing account fails minutes into a queued job.

ISMIP7_SITE_NAME="template"
# Hostname patterns (bash globs, space separated) that identify this cluster.
ISMIP7_SITE_MATCH=""

# --- REQUIRED: software ------------------------------------------------
# The activate script of the Firedrake venv, and the module loads it was
# BUILT against. They must match the build exactly or MPI will mismatch.
ISMIP7_FIREDRAKE="${ISMIP7_FIREDRAKE:-}"
ISMIP7_MODULES="${ISMIP7_MODULES:-}"

# --- REQUIRED: scheduler -----------------------------------------------
# Partitions by job length, the account to charge (empty = your default),
# and a node feature to pin (empty = any node).
ISMIP7_PART_LONG="${ISMIP7_PART_LONG:-}"    # multi-day: 2 km inversions
ISMIP7_PART_SHORT="${ISMIP7_PART_SHORT:-}"  # up to a day: forwards, chained projections
ISMIP7_PART_DEBUG="${ISMIP7_PART_DEBUG:-}"  # minutes: smoke tests and probes
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-}"
ISMIP7_CONSTRAINT="${ISMIP7_CONSTRAINT:-}"

# --- REQUIRED: paths ---------------------------------------------------
# The checkout, and the filesystem that can hold the forcing tree (about 25 GB
# per ESM and scenario, 313 GB for the whole AIS tree). Home directories
# usually cannot; a projects or scratch share can. ISMIP7_REPO_SELF is the
# checkout the submission came from, which is nearly always the right one.
ISMIP7_REPO="${ISMIP7_REPO:-$ISMIP7_REPO_SELF}"
ISMIP7_WORK="${ISMIP7_WORK:-}"

# --- optional: default job size ----------------------------------------
# Ranks per node on PHYSICAL cores (the solver is memory-bandwidth bound, so
# hyperthreads halve per-rank bandwidth), memory per node, wall limits.
ISMIP7_TASKS="${ISMIP7_TASKS:-16}"
ISMIP7_MEM="${ISMIP7_MEM:-120G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-2-00:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-1-00:00:00}"

# Inversions and forwards have different measured sizes. Leave these unset and
# both kinds use the single pair above; set them once this cluster's own
# numbers are known.
#   ISMIP7_TASKS_INV / ISMIP7_MEM_INV    submit.sh inversion
#   ISMIP7_TASKS_FWD / ISMIP7_MEM_FWD    submit.sh projection
#
# ISMIP7_CONSTRAINT splits the same way when the inversions need a different
# node generation from everything else:
#   ISMIP7_CONSTRAINT_INV                submit.sh inversion
#   ISMIP7_CONSTRAINT_FWD                every other kind
