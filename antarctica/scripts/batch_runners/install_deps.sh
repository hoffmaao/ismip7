#!/bin/bash
# Finish the Firedrake venv: the data stack and the three editable ice
# packages. Run ONCE on a login node after build_firedrake.sbatch:
#
#   bash antarctica/scripts/batch_runners/install_deps.sh
#
# The venv, the modules and the work filesystem come from this cluster's
# sites/<name>.sh, the same as every other runner. Set FD_PREFIX if the build
# went somewhere other than $ISMIP7_WORK/sw.
#
# The four sources are rsynced into $FD_PREFIX/src. A clone would miss both of
# the reasons they live there: icepack2 carries two uncommitted functional
# edits (model/minimization.py, model/variational.py) that the inversion
# depends on, and icepack_tools (the shared level-set code) is a local project
# with no remote.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
. "$here/site_env.sh"
ismip7_site_require

PREFIX="${FD_PREFIX:-$ISMIP7_WORK/sw}"
SRC="$PREFIX/src"
export TMPDIR="${TMPDIR:-$ISMIP7_WORK/tmp}"; mkdir -p "$TMPDIR"
ismip7_activate
echo "venv: $VIRTUAL_ENV   python: $(python -V 2>&1)"
echo "src:  $SRC"
python -m pip install --quiet --upgrade pip
python -m pip install --quiet rasterio xarray netCDF4 scipy matplotlib shapely geopandas gmsh pyproj pandas tqdm
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
