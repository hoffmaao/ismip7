# Rice NOTS. Every number here was measured on the login node in September
# 2026; antarctica/scripts/batch_runners/readme.md says how.

ISMIP7_SITE_NAME="rice_nots"
ISMIP7_SITE_MATCH="nots* login*.nots.rice.edu bb*"

# Built by build_firedrake_rice.sbatch (job 1288484). libGLU is in the list because
# the gmsh wheel dlopens libGLU.so.1 at import and the compute nodes lack it.
ISMIP7_FIREDRAKE="${ISMIP7_FIREDRAKE:-/projects/ah301/sw/venv-firedrake/bin/activate}"
ISMIP7_MODULES="${ISMIP7_MODULES:-foss/2023b Python/3.11.5 CMake/3.27.6 M4/1.4.19 flex/2.6.4 Bison/3.8.2 libevent/2.1.12 zstd/1.5.5 Szip/2.1.1 libaec/1.0.6 libGLU/9.0.3}"

# The default account (commons) is granted and needs no -A. It gives commons
# (1 day, 142 nodes), long (3 days, 61 nodes, the only one that fits a 2 km
# inversion), debug (30 min) and scavenge (1 hour, preemptible). The deepsC
# condo (52 uniform Cascade Lake nodes, 187 GB each, no wall limit) needs the
# separate `deepsc` account from the EEPS condo owner, which we do not hold.
ISMIP7_PART_LONG="${ISMIP7_PART_LONG:-long}"
ISMIP7_PART_SHORT="${ISMIP7_PART_SHORT:-commons}"
ISMIP7_PART_DEBUG="${ISMIP7_PART_DEBUG:-scavenge}"
ISMIP7_ACCOUNT="${ISMIP7_ACCOUNT:-}"
# Sapphire Rapids (192 threads, 257-515 GB) is what the 2 km inversions need,
# because `long` holds exactly one Cascade Lake node. Everything else stays on
# cascadelake (80 threads, 187 GB): the measured timings and the chain depth
# come from that generation, and build_firedrake_rice.sbatch configures with
# -march=native, so a build has to land where the runs do.
ISMIP7_CONSTRAINT_INV="${ISMIP7_CONSTRAINT_INV:-sapphirerapids}"
ISMIP7_CONSTRAINT_FWD="${ISMIP7_CONSTRAINT_FWD:-cascadelake}"

# /projects is a 20 TB share; /home is a 10 TB NFS export and /scratch is
# purged, so the forcing tree lives on /projects.
ISMIP7_REPO="${ISMIP7_REPO:-/projects/ah301/ismip7}"
ISMIP7_WORK="${ISMIP7_WORK:-/projects/ah301}"

ISMIP7_TASKS="${ISMIP7_TASKS:-32}"
ISMIP7_MEM="${ISMIP7_MEM:-240G}"
# The forward was measured at 12 ranks and 80 GB resident (readme.md's timing
# table), and the chain depth quoted there is a 12-rank figure. The inversion
# keeps the 32 / 240G above, which is what the 2 km runs need.
ISMIP7_TASKS_FWD="${ISMIP7_TASKS_FWD:-12}"
ISMIP7_MEM_FWD="${ISMIP7_MEM_FWD:-96G}"
ISMIP7_TIME_INV="${ISMIP7_TIME_INV:-3-00:00:00}"
ISMIP7_TIME_FWD="${ISMIP7_TIME_FWD:-1-00:00:00}"
