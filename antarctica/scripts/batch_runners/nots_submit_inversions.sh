#!/bin/bash
# The two 2 km inversion strategies, as a one-factor-at-a-time pair against the
# vertex-sampled 2 km / 20 km-interior run already going on the workstation:
#   B  sampling   : 2 km / 20 km interior, cell_mean   (same mesh as the local run)
#   C  resolution : 2 km /  5 km interior, vertex      (the ~1 M-vertex mesh)
# Both go to `long` on Sapphire Rapids (96 physical cores, 257-515 GB), 32
# ranks. Measured on the workstation, the 20 km-interior mesh sits at 46 GB
# resident at 12 ranks, so 200 G / 240 G are generous, not tight.
#   bash antarctica/scripts/batch_runners/nots_submit_inversions.sh [B|C|BC]
set -euo pipefail
cd "$(dirname "$0")/../../.."
REPO=$PWD
which=${1:-BC}
FRICTION="${ISMIP7_FRICTION:-regularized_coulomb}"
case "$FRICTION" in
  regularized_coulomb) FTAG=_rc ;;
  budd)                FTAG=_budd ;;
  *)                   FTAG= ;;
esac
COMMON="ISMIP7_LC=2000,ISMIP7_BUFFER_M=0,ISMIP7_DHDT_NET_SIGMA=10,ISMIP7_LOG_VEL_WEIGHT=auto,ISMIP7_FRICTION=$FRICTION"
if [[ "$which" == *C* ]]; then
  sbatch -p long -C sapphirerapids --mem=240G -N1 -n32 --cpus-per-task=1 --time=3-00:00:00 -J inv2k_int5k \
    --export=ALL,$COMMON,ISMIP7_LC_COARSE=5000,ISMIP7_MESH=$REPO/antarctica/mesh/antarctica_5000_2000_buffered0.msh,ISMIP7_RASTER_SAMPLE=vertex,ISMIP7_MAP_OUT=$REPO/antarctica/mesh/inversion_icepack2${FTAG}_n3_dg0_logvelnet_2000_int5000.h5 \
    antarctica/scripts/batch_runners/nots_inversion.sbatch
fi
if [[ "$which" == *B* ]]; then
  sbatch -p long -C sapphirerapids --mem=200G -N1 -n32 --cpus-per-task=1 --time=3-00:00:00 -J inv2k_cellmean \
    --export=ALL,$COMMON,ISMIP7_LC_COARSE=20000,ISMIP7_MESH=$REPO/antarctica/mesh/antarctica_20000_2000_buffered0.msh,ISMIP7_RASTER_SAMPLE=cell_mean,ISMIP7_MAP_OUT=$REPO/antarctica/mesh/inversion_icepack2${FTAG}_n3_dg0_logvelnet_2000_cellmean.h5 \
    antarctica/scripts/batch_runners/nots_inversion.sbatch
fi
