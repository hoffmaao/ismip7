#!/bin/bash
# The two 2 km inversion strategies, as a one-factor-at-a-time pair against the
# vertex-sampled 2 km / 20 km-interior run already going on the workstation:
#   B  sampling   : 2 km / 20 km interior, cell_mean   (same mesh as the local run)
#   C  resolution : 2 km /  5 km interior, vertex      (the ~1 M-vertex mesh)
# Both take this cluster's inversion defaults from sites/<cluster>.sh through
# submit.sh (at Rice: long, Sapphire Rapids, 32 ranks). Measured on the workstation, the 20 km-interior mesh sits at 46 GB
# resident at 12 ranks, so 200 G / 240 G are generous, not tight.
#   bash antarctica/scripts/batch_runners/submit_inversions.sh [B|C|BC]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here/../../.."
REPO=$PWD
which=${1:-BC}
. antarctica/scripts/ismip7_names.sh
FRICTION="${ISMIP7_FRICTION:-regularized_coulomb}"
MESH_DIR="$REPO/antarctica/mesh"
COMMON=(ISMIP7_LC=2000 ISMIP7_BUFFER_M=0 ISMIP7_DHDT_NET_SIGMA=10
        ISMIP7_LOG_VEL_WEIGHT=auto "ISMIP7_FRICTION=$FRICTION")
if [[ "$which" == *C* ]]; then
  "$here/submit.sh" inversion --name inv2k_int5k "${COMMON[@]}" \
    ISMIP7_LC_COARSE=5000 \
    "ISMIP7_MESH=$REPO/antarctica/mesh/antarctica_5000_2000_buffered0.msh" \
    ISMIP7_RASTER_SAMPLE=vertex \
    "ISMIP7_MAP_OUT=$MESH_DIR/$(ismip7_map_basename "$FRICTION" 2000_int5000)"
fi
if [[ "$which" == *B* ]]; then
  "$here/submit.sh" inversion --name inv2k_cellmean --mem 200G "${COMMON[@]}" \
    ISMIP7_LC_COARSE=20000 \
    "ISMIP7_MESH=$REPO/antarctica/mesh/antarctica_20000_2000_buffered0.msh" \
    ISMIP7_RASTER_SAMPLE=cell_mean \
    "ISMIP7_MAP_OUT=$MESH_DIR/$(ismip7_map_basename "$FRICTION" 2000_cellmean)"
fi
