#!/bin/bash
# Finish the Firedrake venv: the data stack and the four editable ice
# packages. Run ONCE on a login node after the Firedrake build:
#
#   bash antarctica/scripts/batch_runners/install_deps.sh
#
# The venv, the modules and the work filesystem come from this cluster's
# sites/<name>.sh, the same as every other runner. Set FD_PREFIX if the build
# went somewhere other than $ISMIP7_WORK/sw. TMPDIR is forced under
# ISMIP7_WORK because /tmp is unwritable on many login nodes, which pip only
# discovers part way through an install.
#
# icepack, tlm_adjoint and icepack_tools are rsynced into $FD_PREFIX/src;
# icepack_tools (the shared level-set code) is a local project with no remote,
# which is why that is the route. icepack2 is NOT rsynced: it is cloned at the
# pin below, so that every site installs one known tree. Run with
# --check-icepack2 to verify that pin and do nothing else.
set -euo pipefail

# --- the icepack2 the inversion needs -----------------------------------
# The inversion needs two functional lines that Firedrake 2026 forces
# (Mesh.geometric_dimension became an attribute, and viscous_power and flow_law
# call it). They lived as uncommitted edits in one workstation checkout, which
# is why the source used to be rsynced and why no reader could tell whether two
# sites ran the same icepack2 (issue #46). They are now a commit on a branch:
# icepack/icepack2#3, open at the time of writing. The pin is that branch's head
# by SHA, which is icepack2 main (40e848b) plus those two files, so a clone
# gives a site exactly what the workstations ran.
#
# When the pull request merges, follow it by overriding all three:
#   ICEPACK2_REMOTE=https://github.com/icepack/icepack2.git ICEPACK2_REF=main \
#   ICEPACK2_SHA=<the commit> bash install_deps.sh --check-icepack2
# and change the defaults here in the same pass.
ICEPACK2_REMOTE="${ICEPACK2_REMOTE:-https://github.com/hoffmaao/icepack2.git}"
ICEPACK2_REF="${ICEPACK2_REF:-fix/firedrake-2026-geometric-dimension}"
ICEPACK2_SHA="${ICEPACK2_SHA:-e0a46c9ce3e95dc916660a200e695f0417aa3037}"

# Clone the pin, or verify the checkout already there is it. A tree that is not
# a clone, or that carries local changes, is reported and left alone: it may be
# the rsynced copy whose edits this pin replaces, and discarding those silently
# would throw away the only copy of whatever else is in them.
ismip7_icepack2_pin() {
    local dir="$SRC/icepack2" head
    if [ -d "$dir" ] && [ ! -d "$dir/.git" ]; then
        echo "ERROR: $dir is not a git checkout. It is the rsynced copy this" >&2
        echo "       pin replaces: move it aside (mv '$dir' '$dir.rsynced')" >&2
        echo "       and run this script again to clone $ICEPACK2_SHA." >&2
        return 2
    fi
    command -v git >/dev/null 2>&1 || {
        echo "ERROR: no git here, so the icepack2 pin cannot be checked." >&2
        echo "       Run this on a login node; compute nodes have no git." >&2
        return 2
    }
    if [ ! -d "$dir" ]; then
        mkdir -p "$SRC"
        echo "icepack2: cloning $ICEPACK2_REMOTE"
        git clone --quiet "$ICEPACK2_REMOTE" "$dir" || return 2
    fi
    if [ -n "$(git -C "$dir" status --porcelain)" ]; then
        echo "ERROR: $dir has local changes, so it is not the pinned tree." >&2
        echo "       Commit or discard them, or move the checkout aside, then" >&2
        echo "       run this script again." >&2
        return 2
    fi
    head="$(git -C "$dir" rev-parse HEAD)"
    if [ "$head" != "$ICEPACK2_SHA" ]; then
        git -C "$dir" fetch --quiet "$ICEPACK2_REMOTE" "$ICEPACK2_REF" \
            || git -C "$dir" fetch --quiet --all || return 2
        git -C "$dir" checkout --quiet --detach "$ICEPACK2_SHA" || {
            echo "ERROR: $ICEPACK2_SHA is not in $dir after fetching" >&2
            echo "       $ICEPACK2_REF from $ICEPACK2_REMOTE." >&2
            return 2
        }
    fi
    echo "icepack2: $ICEPACK2_SHA ($ICEPACK2_REF from $ICEPACK2_REMOTE)"
}
here="$(cd "$(dirname "$0")" && pwd)"
. "$here/site_env.sh"
ismip7_site_require

PREFIX="${FD_PREFIX:-$ISMIP7_WORK/sw}"
SRC="$PREFIX/src"

# The pin alone, for a site checking what it runs without touching its venv.
if [ "${1:-}" = "--check-icepack2" ]; then
    ismip7_icepack2_pin || exit $?
    exit 0
fi

export TMPDIR="$ISMIP7_WORK/tmp"; mkdir -p "$TMPDIR"
ismip7_activate
echo "venv: $VIRTUAL_ENV   python: $(python -V 2>&1)"
echo "src:  $SRC"
python -m pip install --quiet --upgrade pip
python -m pip install --quiet rasterio xarray netCDF4 scipy matplotlib shapely geopandas gmsh pyproj pandas tqdm
ismip7_icepack2_pin || exit $?
for d in icepack icepack2 tlm_adjoint icepack_tools; do
    [ -d "$SRC/$d" ] || { echo "ERROR: $SRC/$d missing (rsync it)"; exit 2; }
    python -m pip install --quiet -e "$SRC/$d"
done
python - <<'PY'
import h5py, netCDF4
print("hdf5: h5py", h5py.version.hdf5_version, "| netCDF4", netCDF4.__hdf5libversion__)
import icepack, icepack2, tlm_adjoint, icepack_tools, rasterio, xarray, gmsh, geopandas, shapely, scipy
print("deps OK: icepack", icepack.__version__)
PY
