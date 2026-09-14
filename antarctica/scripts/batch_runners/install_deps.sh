#!/bin/bash
# Finish the NOTS Firedrake venv: the data stack and the three editable ice
# packages. Run ONCE on a login node after build_firedrake.sbatch:
#
#   bash antarctica/scripts/batch_runners/install_deps.sh
#
# Sources are rsynced from the workstation into /projects/ah301/sw/src, not
# cloned, because icepack2 carries two uncommitted functional edits
# (model/minimization.py, model/variational.py) the inversion depends on, and
# icepack_tools (the shared level-set code) is a local project with no remote.
set -euo pipefail
export TMPDIR=/projects/ah301/tmp; mkdir -p "$TMPDIR"     # /tmp is not writable here
. /projects/ah301/ismip7/antarctica/scripts/batch_runners/site_env.sh
ismip7_activate
echo "venv: $VIRTUAL_ENV   python: $(python -V 2>&1)"
python -m pip install --quiet --upgrade pip
python -m pip install --quiet rasterio xarray netCDF4 scipy matplotlib shapely geopandas gmsh pyproj pandas tqdm
for d in icepack icepack2 tlm_adjoint icepack_tools; do
    [ -d "/projects/ah301/sw/src/$d" ] || { echo "ERROR: /projects/ah301/sw/src/$d missing (rsync it)"; exit 2; }
    python -m pip install --quiet -e "/projects/ah301/sw/src/$d"
done
python - <<'PY'
import h5py, netCDF4
print("hdf5: h5py", h5py.version.hdf5_version, "| netCDF4", netCDF4.__hdf5libversion__)
import icepack, icepack2, tlm_adjoint, icepack_tools, rasterio, xarray, gmsh, geopandas, shapely, scipy
print("deps OK: icepack", icepack.__version__)
PY
