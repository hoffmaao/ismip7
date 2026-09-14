# University of Chicago (RCC Midway). NOT yet filled in: nobody has run this
# pipeline there, so rather than guess partition names, an account and a
# module stack, the required fields are left empty and site_env.sh refuses to
# submit until they are set. Fill them from RCC's own documentation and from
# `sinfo -s` / `accounts balance` on the login node, then delete this notice.
#
# What each field needs:
#   ISMIP7_FIREDRAKE  the activate script of a Firedrake 2026.4 build there
#                     (build_firedrake.sbatch is the Rice recipe; the module
#                     names will differ, the sequence will not)
#   ISMIP7_MODULES    the loads that build was made against
#   ISMIP7_PART_*     partitions by job length (`sinfo -s`)
#   ISMIP7_ACCOUNT    the allocation to charge, usually pi-<name>
#   ISMIP7_REPO       the checkout; ISMIP7_WORK a share that can hold about
#                     25 GB of forcing per ESM and scenario, 313 GB for the
#                     whole AIS tree (home usually cannot)

ISMIP7_SITE_NAME="uchicago_midway"
ISMIP7_SITE_MATCH="midway* login*.rcc.uchicago.edu"

ISMIP7_FIREDRAKE="${ISMIP7_FIREDRAKE:-}"
ISMIP7_MODULES="${ISMIP7_MODULES:-}"

ISMIP7_PART_LONG="${ISMIP7_PART_LONG:-}"
ISMIP7_PART_SHORT="${ISMIP7_PART_SHORT:-}"
ISMIP7_PART_DEBUG="${ISMIP7_PART_DEBUG:-}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-}"
ISMIP7_CONSTRAINT="${ISMIP7_CONSTRAINT:-}"

ISMIP7_REPO="${ISMIP7_REPO:-$HOME/ismip7}"
ISMIP7_WORK="${ISMIP7_WORK:-}"

ISMIP7_TASKS="${ISMIP7_TASKS:-16}"
ISMIP7_MEM="${ISMIP7_MEM:-120G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-2-00:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-1-00:00:00}"
