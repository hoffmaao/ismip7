# University of Chicago RCC Midway3. Firedrake here is a container image, not a
# venv: icepack2_midway3_source.def at the repository root builds it
# (Firedrake 2026.4.1, PETSc 3.25.5 with parmetis). Set ISMIP7_CONTAINER to the
# .sif and ISMIP7_FIREDRAKE is not asked for; site_core.sh then runs every
# `python` in the image and starts MPI with the image's own mpiexec.
#
# NOT yet run end to end: nobody has submitted this pipeline at RCC. Rather
# than guess partition names and an account, those fields are empty and
# submit.sh refuses until they are set, here or in sites/local.env. What each
# needs, and where to read it on a login node (site_recon.sh prints most):
#   ISMIP7_CONTAINER   the built image, e.g. /project/<pi>/ismip7/icepack2.sif
#                      (per user or per group: sites/local.env is the place)
#   ISMIP7_MODULES     whatever puts `apptainer` on PATH on a compute node, if
#                      anything is needed (`module avail apptainer singularity`)
#   ISMIP7_PART_*      partitions by job length (`sinfo -s`)
#   ISMIP7_ACCOUNT     the allocation to charge, usually pi-<name>
#                      (`accounts balance`)
#   ISMIP7_WORK        a share that can hold about 25 GB of forcing per ESM and
#                      scenario, 313 GB for the whole AIS tree (home cannot)
#   ISMIP7_CORES_PER_NODE / ISMIP7_MEM_PER_NODE
#                      one node of the partition the timing lanes will use
#                      (`sinfo -p <partition> -N -o "%c %m %f" | sort -u`);
#                      lanes that do not fit are recorded not runnable
# Then: submit.sh verify (4 ranks; checks the image sees the checkout and that
# MPI works across ranks in it), submit.sh smoke, make -C antarctica timing-dry-run.
# If mpiexec in the image objects to running under a Slurm allocation, try
#   ISMIP7_CONTAINER_MPIEXEC="mpiexec --mca plm isolated"
# and delete this notice once verify passes.

ISMIP7_SITE_NAME="uchicago_midway"
ISMIP7_SITE_MATCH="midway* login*.rcc.uchicago.edu"

ISMIP7_CONTAINER="${ISMIP7_CONTAINER:-}"
ISMIP7_CONTAINER_RUNTIME="${ISMIP7_CONTAINER_RUNTIME:-apptainer}"
ISMIP7_CONTAINER_ARGS="${ISMIP7_CONTAINER_ARGS:-}"
ISMIP7_FIREDRAKE="${ISMIP7_FIREDRAKE:-}"
ISMIP7_MODULES="${ISMIP7_MODULES:-}"

ISMIP7_PART_LONG="${ISMIP7_PART_LONG:-}"
ISMIP7_PART_SHORT="${ISMIP7_PART_SHORT:-}"
ISMIP7_PART_DEBUG="${ISMIP7_PART_DEBUG:-}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-}"
ISMIP7_CONSTRAINT="${ISMIP7_CONSTRAINT:-}"

ISMIP7_REPO="${ISMIP7_REPO:-$ISMIP7_REPO_SELF}"
ISMIP7_WORK="${ISMIP7_WORK:-}"

# The image keeps its kernel caches under /scratch/midway3/$USER only when it
# is started with `apptainer run`; the runners use `exec`, so the timing
# lanes' persistent cache is named here instead.
ISMIP7_TIMING_JIT_CACHE="${ISMIP7_TIMING_JIT_CACHE:-/scratch/midway3/${USER:-$(id -un)}/icepack_cache}"

ISMIP7_TASKS="${ISMIP7_TASKS:-16}"
ISMIP7_MEM="${ISMIP7_MEM:-120G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-2-00:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-1-00:00:00}"

ISMIP7_CORES_PER_NODE="${ISMIP7_CORES_PER_NODE:-}"
ISMIP7_MEM_PER_NODE="${ISMIP7_MEM_PER_NODE:-}"
