# "Whatever is already in the environment." Never matched by hostname: ask for
# it with ISMIP7_SITE=local. This is how the runner logic is exercised off a
# cluster (tests/test_inversion_chain.py, tests/test_projection_chain.py) and
# how you debug a job script on a workstation, where the venv is already
# active and there is no scheduler to satisfy.

ISMIP7_SITE_NAME="local"
ISMIP7_SITE_MATCH=""

ISMIP7_FIREDRAKE="${ISMIP7_FIREDRAKE:-${VIRTUAL_ENV:+$VIRTUAL_ENV/bin/activate}}"
ISMIP7_MODULES="${ISMIP7_MODULES:-}"

ISMIP7_PART_LONG="${ISMIP7_PART_LONG:-local}"
ISMIP7_PART_SHORT="${ISMIP7_PART_SHORT:-local}"
ISMIP7_PART_DEBUG="${ISMIP7_PART_DEBUG:-local}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-}"
ISMIP7_CONSTRAINT="${ISMIP7_CONSTRAINT:-}"

ISMIP7_REPO="${ISMIP7_REPO:-$PWD}"
ISMIP7_WORK="${ISMIP7_WORK:-$ISMIP7_REPO}"

ISMIP7_TASKS="${ISMIP7_TASKS:-4}"
ISMIP7_MEM="${ISMIP7_MEM:-16G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-01:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-01:00:00}"
