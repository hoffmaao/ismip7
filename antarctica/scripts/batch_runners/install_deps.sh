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
# icepack and icepack_tools are rsynced into $FD_PREFIX/src; icepack_tools
# (the shared level-set code) is a local project with no remote, which is why
# that is the route. icepack2 and tlm_adjoint are NOT rsynced: each is cloned
# at the pin below, so that every site installs one known tree. Run with
# --check-icepack2 or --check-tlm-adjoint to verify that pin and do nothing
# else.
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

# --- the tlm_adjoint the inversion and the uncertainty runs need ---------
# Upstream tlm_adjoint (0.9.1dev, April 2026) leaves each Gauss-Newton
# action's tangent-linear copies of the forward variables attached to them:
# one copy of the forward state per action, 0.6 GB per Hessian action at 2 km,
# which ended a 100-mode eigensolve out of memory. Our fork carries the fix on
# its ismip7 branch, pinned here by SHA like icepack2; tests/test_laplace.py
# fails on an environment that still installs upstream. Follow the fork the
# same way, overriding all three and changing the defaults in the same pass.
TLM_ADJOINT_REMOTE="${TLM_ADJOINT_REMOTE:-https://github.com/hoffmaao/tlm_adjoint.git}"
TLM_ADJOINT_REF="${TLM_ADJOINT_REF:-ismip7}"
TLM_ADJOINT_SHA="${TLM_ADJOINT_SHA:-e9c970e2f0786119aa4d41410b5dbf54ac546b75}"

# Clone a pin into $SRC/<name>, or verify the checkout already there is it. A
# tree that is not a clone, or that carries local changes to tracked files, is
# reported and left alone: it may be the rsynced copy whose edits this pin
# replaces, and discarding those silently would throw away the only copy of
# whatever else is in them. Untracked files (an editable install's egg-info)
# are not local changes.
ismip7_git_pin() {
    local name="$1" remote="$2" ref="$3" sha="$4"
    local dir="$SRC/$name" head
    if [ -d "$dir" ] && [ ! -d "$dir/.git" ]; then
        echo "ERROR: $dir is not a git checkout. It is the rsynced copy this" >&2
        echo "       pin replaces: move it aside (mv '$dir' '$dir.rsynced')" >&2
        echo "       and run this script again to clone $sha." >&2
        return 2
    fi
    command -v git >/dev/null 2>&1 || {
        echo "ERROR: no git here, so the $name pin cannot be checked." >&2
        echo "       Run this on a login node; compute nodes have no git." >&2
        return 2
    }
    if [ ! -d "$dir" ]; then
        mkdir -p "$SRC"
        echo "$name: cloning $remote"
        git clone --quiet "$remote" "$dir" || return 2
    fi
    if [ -n "$(git -C "$dir" status --porcelain --untracked-files=no)" ]; then
        echo "ERROR: $dir has local changes, so it is not the pinned tree." >&2
        echo "       Commit or discard them, or move the checkout aside, then" >&2
        echo "       run this script again." >&2
        return 2
    fi
    head="$(git -C "$dir" rev-parse HEAD)"
    if [ "$head" != "$sha" ]; then
        git -C "$dir" fetch --quiet "$remote" "$ref" \
            || git -C "$dir" fetch --quiet --all || return 2
        git -C "$dir" checkout --quiet --detach "$sha" || {
            echo "ERROR: $sha is not in $dir after fetching" >&2
            echo "       $ref from $remote." >&2
            return 2
        }
    fi
    echo "$name: $sha ($ref from $remote)"
}
ismip7_icepack2_pin() {
    ismip7_git_pin icepack2 "$ICEPACK2_REMOTE" "$ICEPACK2_REF" "$ICEPACK2_SHA"
}
ismip7_tlm_adjoint_pin() {
    ismip7_git_pin tlm_adjoint "$TLM_ADJOINT_REMOTE" "$TLM_ADJOINT_REF" "$TLM_ADJOINT_SHA"
}
here="$(cd "$(dirname "$0")" && pwd)"
. "$here/site_env.sh"
ismip7_site_require

PREFIX="${FD_PREFIX:-$ISMIP7_WORK/sw}"
SRC="$PREFIX/src"

# A pin alone, for a site checking what it runs without touching its venv.
case "${1:-}" in
    --check-icepack2) ismip7_icepack2_pin || exit $?; exit 0 ;;
    --check-tlm-adjoint) ismip7_tlm_adjoint_pin || exit $?; exit 0 ;;
esac

export TMPDIR="$ISMIP7_WORK/tmp"; mkdir -p "$TMPDIR"
ismip7_activate
echo "venv: $VIRTUAL_ENV   python: $(python -V 2>&1)"
echo "src:  $SRC"
python -m pip install --quiet --upgrade pip
python -m pip install --quiet rasterio xarray netCDF4 scipy matplotlib shapely geopandas gmsh pyproj pandas tqdm
ismip7_icepack2_pin || exit $?
ismip7_tlm_adjoint_pin || exit $?
for d in icepack icepack_tools; do
    [ -d "$SRC/$d" ] || { echo "ERROR: $SRC/$d missing (rsync it)"; exit 2; }
done
for d in icepack icepack2 tlm_adjoint icepack_tools; do
    python -m pip install --quiet -e "$SRC/$d"
done
python - <<'PY'
import h5py, netCDF4
print("hdf5: h5py", h5py.version.hdf5_version, "| netCDF4", netCDF4.__hdf5libversion__)
import icepack, icepack2, tlm_adjoint, icepack_tools, icepack_tools.calving, rasterio, xarray, gmsh, geopandas, shapely, scipy
print("deps OK: icepack", icepack.__version__)
PY
