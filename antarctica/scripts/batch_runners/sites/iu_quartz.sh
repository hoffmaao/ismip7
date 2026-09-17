# Indiana University Quartz. Taken from the IU runners on the upstream
# `timing_matrix` branch (antarctica/scripts/batch_runners/timing_*.script),
# which is where these module names, the account and the venv path come from.
# The defaults below point at an existing IU build. A second IU user keeps the
# module stack and, in sites/local.env, points ISMIP7_FIREDRAKE at their own
# build (or at that one, if it is readable to them) and ISMIP7_ACCOUNT at
# their own allocation.

ISMIP7_SITE_NAME="iu_quartz"
ISMIP7_SITE_MATCH="quartz* h2.quartz* login*.quartz.uits.iu.edu"

ISMIP7_FIREDRAKE="${ISMIP7_FIREDRAKE:-/N/u/dlilien/Quartz/sw/firedrake/2026.04/firedrake/bin/activate}"
# `module use` first: the firedrake and petsc modulefiles are in a personal
# tree, not the system one.
ISMIP7_MODULES="${ISMIP7_MODULES:-gnu/9.3.0 openmpi/4.0.5 python/3.14.5 zlib/1.2.13 hdf5/1.12 openblas/0.3.13 patchelf/0.18.0 petsc/3.25.5 firedrake/2026.04}"
ISMIP7_MODULE_USE="${ISMIP7_MODULE_USE:-/N/u/dlilien/Quartz/modulefiles}"

ISMIP7_PART_LONG="${ISMIP7_PART_LONG:-general}"
ISMIP7_PART_SHORT="${ISMIP7_PART_SHORT:-general}"
ISMIP7_PART_DEBUG="${ISMIP7_PART_DEBUG:-debug}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-r00905}"
ISMIP7_CONSTRAINT="${ISMIP7_CONSTRAINT:-}"

# The checkout the submission came from (site_core.sh works it out), so a
# second worktree runs its own tree without anything being edited.
ISMIP7_REPO="${ISMIP7_REPO:-$ISMIP7_REPO_SELF}"
ISMIP7_WORK="${ISMIP7_WORK:-$HOME}"

# The IU timing runs use 12 to 16 ranks per node, 128 to 240 GB, up to 48 h.
ISMIP7_TASKS="${ISMIP7_TASKS:-16}"
ISMIP7_MEM="${ISMIP7_MEM:-240G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-48:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-24:00:00}"
