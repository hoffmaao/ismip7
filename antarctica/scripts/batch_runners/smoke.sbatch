#!/bin/bash
# Smoke test: the 32 km inversion for 2 iterates on 4 ranks. Exercises the
# whole stack (firedrake, icepack, icepack2, tlm_adjoint, rasterio, every data
# file) in a few minutes, so a multi-day job cannot die on a missing file or
# module hours after it queued. It is also the first thing to run at a new
# site, once sites/<name>.sh is filled in.
#   antarctica/scripts/batch_runners/submit.sh smoke
# Resources (partition, tasks, memory, wall time, account, constraint) come
# from submit.sh and this cluster's sites/<name>.sh, never from directives
# here: a directive cannot read a site file, and a header that disagrees
# with the command line is rejected outright by Slurm.
#SBATCH -o logs/ismip7_smoke_%j.out
#SBATCH -e logs/ismip7_smoke_%j.err
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
# Defaults are the 32 km smoke; every value can be overridden with a KEY=VALUE
# on the submit.sh line to smoke-test another mesh or sampling, e.g. two
# iterates of the 2 km cell-mean case on 32 Sapphire Rapids ranks. They MUST
# be set before site_env.sh is sourced: it fills ISMIP7_LC / ISMIP7_MESH with
# the 2500 m production defaults, after which a ${VAR:-...} fallback here never
# fires (job 1300996 silently ran the 2500 m mesh on 4 ranks that way).
export ISMIP7_LC="${ISMIP7_LC:-32000}" ISMIP7_LC_COARSE="${ISMIP7_LC_COARSE:-320000}" ISMIP7_BUFFER_M="${ISMIP7_BUFFER_M:-0}"
# The mesh is not tracked, so a new site has none. Without ISMIP7_MESH the smoke
# takes the mesh under the name the generator writes today, or the unsuffixed
# file an older generator left (what Rice has), and builds the first when
# neither is there: a minute of serial gmsh, once per checkout. A mesh named
# through ISMIP7_MESH is the caller's and is never built here.
smoke_mesh_dir="${ISMIP7_REPO:-.}/antarctica/mesh"
smoke_mesh="$smoke_mesh_dir/antarctica_${ISMIP7_LC_COARSE}_${ISMIP7_LC}_buffered${ISMIP7_BUFFER_M%.*}.msh"
smoke_mesh_old="$smoke_mesh_dir/antarctica_${ISMIP7_LC_COARSE}_${ISMIP7_LC}.msh"
smoke_build_mesh=""
if [ -z "${ISMIP7_MESH:-}" ]; then
    if [ ! -f "$smoke_mesh" ] && [ -f "$smoke_mesh_old" ]; then
        smoke_mesh="$smoke_mesh_old"
    elif [ ! -f "$smoke_mesh" ]; then
        smoke_build_mesh=1
    fi
    export ISMIP7_MESH="$smoke_mesh"
fi
. antarctica/scripts/batch_runners/site_env.sh
ismip7_activate
ismip7_banner
if [ -n "$smoke_build_mesh" ]; then
    echo "=== $(date -Is) no $ISMIP7_MESH; building it ==="
    python -u antarctica/scripts/mesh_antarctica.py --lc "$ISMIP7_LC" \
        --lc-coarse "$ISMIP7_LC_COARSE" --buffer-m "$ISMIP7_BUFFER_M" \
        || { rc=$?; echo "=== $(date -Is) smoke mesh build exit $rc ==="; exit $rc; }
fi
export ISMIP7_MISFIT_NORM=sigma ISMIP7_SIGMA_U_FLOOR=3 ISMIP7_LOG_VEL_WEIGHT=auto
export ISMIP7_DHDT_WEIGHT=1 ISMIP7_DHDT_NET_SIGMA=10 ISMIP7_GAMMA_THETA=1e5 ISMIP7_GAMMA_PHI=1e5
export ISMIP7_MAXITER="${ISMIP7_MAXITER:-2}" ISMIP7_RASTER_SAMPLE="${ISMIP7_RASTER_SAMPLE:-vertex}"
export ISMIP7_MAP_OUT="${ISMIP7_MAP_OUT:-$ISMIP7_REPO/antarctica/results/smoke_${ISMIP7_LC}_${ISMIP7_LC_COARSE}_${ISMIP7_RASTER_SAMPLE}.h5}"
mkdir -p "$(dirname "$ISMIP7_MAP_OUT")"
ismip7_mpirun "${SLURM_NTASKS:-1}" python -u antarctica/scripts/inversion_icepack2.py
rc=$?
echo "=== $(date -Is) smoke exit $rc ==="
exit $rc
