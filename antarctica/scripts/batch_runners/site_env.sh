#!/bin/bash
# Shared environment for ISMIP7 icepack2 batch jobs, on any cluster.
# Sourced by the inversion and projection runners and by submit.sh. Nothing
# here submits anything.
#
# Two halves. site_core.sh answers the site questions (which cluster, where is
# Firedrake, what does the scheduler want, the per-user sites/local.env) and
# provides ismip7_activate, ismip7_chain_resources and ismip7_banner. This file
# adds the model configuration every inversion and projection shares. A job
# that brings its whole model configuration with it, as the timing campaign's
# lanes do, sources site_core.sh alone.
#
#   ISMIP7_SITE=iu_quartz sbatch ...     # name it
#   sbatch ...                           # or let the hostname choose
# shellcheck disable=SC1091
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/site_core.sh"

# --- repository and data ------------------------------------------------
# Two roots, because a cluster can hold more than one checkout. ISMIP7_REPO is
# the code, and follows the checkout the command was run from. ISMIP7_SHARE is
# the tree holding the large artifacts no checkout carries -- the forcing tree,
# the meshes, the MAPs and the observational rasters, all gitignored -- which a
# site file names when they do not live beside the code. It defaults to
# ISMIP7_REPO, so a site holding one checkout behaves as it always has.
export ISMIP7_SHARE="${ISMIP7_SHARE:-$ISMIP7_REPO}"
export ISMIP7_DATA_ROOT="${ISMIP7_DATA_ROOT:-$ISMIP7_SHARE/ISMIP7/AIS}"
# simulation.py and preflight.py default this to the running checkout's
# antarctica/data, which is gitignored and empty in a fresh clone, so state it
# here too: the shell layer and the Python layer then name the same rasters.
export ISMIP7_OBS_DATA_ROOT="${ISMIP7_OBS_DATA_ROOT:-$ISMIP7_SHARE/antarctica/data}"

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
export ISMIP7_MESH="${ISMIP7_MESH:-$ISMIP7_SHARE/antarctica/mesh/antarctica_64000_2500.msh}"
# The MAP this configuration writes and reads: named from the law, so switching
# ISMIP7_FRICTION switches the file and a Budd re-inversion cannot land on the
# RC MAP. One variable for both halves of the workflow - inversion.sbatch
# writes it, projection.sbatch loads it - because the runners' `logvelnet`
# name is not one a forward can derive for itself.
. "$_ISMIP7_BR_DIR/../ismip7_names.sh"
export ISMIP7_MAP_DEFAULT="${ISMIP7_MAP_DEFAULT:-$ISMIP7_SHARE/antarctica/mesh/$(ismip7_map_basename "$ISMIP7_FRICTION" "$ISMIP7_LC")}"

ismip7_banner_model() {
    echo "    mesh    $ISMIP7_MESH"
    echo "    lc=$ISMIP7_LC lc_coarse=$ISMIP7_LC_COARSE geometry=$ISMIP7_GEOMETRY_SPACE" \
         "friction=$ISMIP7_FRICTION n=$ISMIP7_N_FLOW"
}
