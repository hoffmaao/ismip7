#!/bin/bash
# Build Firedrake + PETSc ONCE, into a shared prefix. This is the RICE recipe:
# the prefix and the module names below are Rice's, so at another site copy it
# and substitute that cluster's own module stack. The sequence does not change.
#
#   antarctica/scripts/batch_runners/submit.sh build
#
# THE POINT OF BUILDING IT HERE RATHER THAN COPYING THE WORKSTATION BUILD:
# the workstation PETSc has NO parallel graph partitioner. Verified 2026-09-12
# by forcing each type through a real 4-rank distribution:
#
#     parmetis  -> PETSc error 56 (not supported: backend absent)
#     ptscotch  -> PETSc error 56
#     chaco     -> PETSc error 56
#     simple    -> works
#
# The workstation's configure line does list `--with-ptscotch`, so the name
# alone proves nothing; the backend was never actually built. PETSc then falls
# back to `simple`, which cuts the DMPlex by point index. On a structured mesh
# that is harmless (ghost/owned 0.11 at 4 ranks). On a gmsh-ordered unstructured
# mesh, index order is spatially arbitrary, and the Antarctica mesh measures
# ghost/owned 16.0 at 4 ranks rising to 287 at 32 -- a halo that DIVERGES with
# rank count, so any scaling number taken on it measures the partitioner rather
# than the model.
#
# WHAT ACTUALLY FIXES IT: upstream's own option set ALREADY carries
# --download-ptscotch (verified 2026-09-12 by reading what firedrake-configure
# emits). So a stock Firedrake build here fixes the partitioner on its own; the
# workstation build is broken because it was configured against a system scotch
# that never provided the parallel partitioner, not because upstream omits it.
# This script adds only --download-parmetis, to have a second partitioner to
# cross-check against. Everything else mirrors the upstream procedure at
# https://www.firedrakeproject.org/install.html.
#
# NOTS is RHEL 9.4, which firedrake-configure does not recognise, so every
# invocation passes --os unknown. That makes PETSc build its dependencies from
# source rather than expecting distro packages, which is what you want inside a
# module environment anyway.
#
# Runs on a COMPUTE node, constrained to cascadelake, deliberately: the configure
# uses -march=native, so the build must happen on the same generation it will run
# on. Building on the login node would bake in whatever that node happens to be.

# Resources (partition, tasks, memory, wall time, account, constraint) come
# from submit.sh and this cluster's sites/<name>.sh, never from directives
# here: a directive cannot read a site file, and a header that disagrees
# with the command line is rejected outright by Slurm.
#SBATCH -o logs/fd_build_%j.out
#SBATCH -e logs/fd_build_%j.err
#SBATCH --mail-type=END,FAIL

set -uo pipefail

PREFIX="${FD_PREFIX:-/projects/ah301/sw}"
PETSC_SRC="$PREFIX/petsc"
VENV="$PREFIX/venv-firedrake"
NPROC="${SLURM_CPUS_PER_TASK:-8}"

echo "=== $(date -Is)  firedrake build on $(hostname) ==="
echo "    prefix   $PREFIX"
echo "    build cpus $NPROC"

# --- modules ------------------------------------------------------------
# foss/2023b is the toolchain that actually resolves here (GCC 13.2.0,
# OpenMPI 4.1.6, OpenBLAS, ScaLAPACK, FFTW). Python 3.11.5 is built against the
# same GCCcore/13.2.0, so the ABI matches. Whatever loads here MUST be recorded
# in site_env.sh as ISMIP7_MODULES, or the runtime will mismatch this build.
module purge
# flex/Bison/M4 are NOT on the compute nodes at all (verified 2026-09-12:
# flex, lex, bison and m4 are all MISSING on bb2u12c1), and PETSc's
# --download-bison refuses to build without flex in PATH. The modules supply
# them. This is the second compute-node gap after git; assume nothing from the
# bare node image.
# Versions pinned explicitly: bare `flex` resolves against the wrong GCCcore
# and silently leaves M4 unloaded, which fails the flex module and drops the
# whole load back to /usr/bin (gcc 11, python 3.9, no mpicc).
# libevent is required even though nothing here asks for it directly: the
# EasyBuild OpenMPI 4.1.6 wrapper puts -levent_core -levent_pthreads on every
# link line, but does not pull libevent onto LIBRARY_PATH. Without it PnetCDF
# (and anything else that link-tests with mpicc) fails with
#   ld: cannot find -levent_core
# which reads like a PnetCDF bug and is not one.
# zstd/Szip/libaec are the same story as libevent: HDF5 and NetCDF detect the
# features from headers the toolchain provides, then link -lzstd / -lsz, but
# nothing puts those libraries on LIBRARY_PATH. NetCDF fails with
#   ld: cannot find -lzstd
# All of these were found by failing one at a time; they are loaded together
# now so the next new package does not cost another 16-minute round trip.
module load foss/2023b Python/3.11.5 CMake/3.27.6 \
            M4/1.4.19 flex/2.6.4 Bison/3.8.2 libevent/2.1.12 \
            zstd/1.5.5 Szip/2.1.1 libaec/1.0.6 || {
    echo "ERROR: module load failed" >&2; exit 2; }
# Fail loudly if the load silently degraded to system tools.
for t in mpicc flex bison; do
    command -v "$t" >/dev/null || { echo "ERROR: $t missing after module load" >&2; exit 2; }
done
case "$(command -v gcc)" in /opt/apps/*) ;; *) echo "ERROR: gcc is $(command -v gcc), modules did not take" >&2; exit 2;; esac
module list 2>&1 | sed 's/^/    /'

echo "    gcc  $(gcc -dumpversion)"
echo "    mpi  $(ompi_info --version 2>/dev/null | head -1)"
echo "    py   $(python3 --version 2>&1)"

mkdir -p "$PREFIX" || { echo "ERROR: cannot create $PREFIX" >&2; exit 2; }
cd "$PREFIX"

# --- firedrake-configure: the authoritative source of PETSc options ------
echo; echo "=== fetching firedrake-configure ==="
curl -fsSL -O https://raw.githubusercontent.com/firedrakeproject/firedrake/release/scripts/firedrake-configure || {
    echo "ERROR: could not fetch firedrake-configure (no egress?)" >&2; exit 3; }

FDC="python3 $PREFIX/firedrake-configure --os unknown"
PETSC_VERSION=$($FDC --show-petsc-version) || exit 3
echo "    petsc version: $PETSC_VERSION"

# --- PETSc ---------------------------------------------------------------
# NO GIT ON THE COMPUTE NODES. Verified 2026-09-12: /usr/bin/git exists on the
# login node but not on bb2u12c1, and no git module is published. So the source
# comes from the GitLab release tarball over curl, which is available (the
# cURL/8.3.0 module rides in with foss/2023b). Do not "fix" this back to a
# git clone; it will fail the moment it leaves the login node.
if [ ! -d "$PETSC_SRC" ]; then
    echo; echo "=== fetching petsc $PETSC_VERSION tarball ==="
    TARBALL="$PREFIX/petsc-${PETSC_VERSION}.tar.gz"
    if [ ! -f "$TARBALL" ]; then
        curl -fL --retry 3 -o "$TARBALL" \
          "https://gitlab.com/petsc/petsc/-/archive/${PETSC_VERSION}/petsc-${PETSC_VERSION}.tar.gz" \
          || { echo "ERROR: petsc tarball download failed" >&2; exit 3; }
    fi
    mkdir -p "$PETSC_SRC"
    tar -xzf "$TARBALL" -C "$PETSC_SRC" --strip-components=1 \
        || { echo "ERROR: petsc tarball extract failed" >&2; exit 3; }
    echo "    extracted to $PETSC_SRC"
fi

export PETSC_DIR="$PETSC_SRC"
export PETSC_ARCH=arch-firedrake-default

if [ ! -f "$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables" ]; then
    echo; echo "=== configuring petsc ==="
    cd "$PETSC_SRC"
    # The option string contains QUOTED flags with embedded spaces, e.g.
    #   --COPTFLAGS='-O3 -march=native -mtune=native'
    # so it must NOT be split on whitespace. `eval` preserves the quoting;
    # splitting it (an earlier version of this script used `tr ' ' '\n'`) would
    # hand configure a broken -O3 and three stray arguments.
    PETSC_OPTS=$($FDC --show-petsc-configure-options -- --download-parmetis) || exit 4
    echo "    configure options:"; echo "$PETSC_OPTS" | fold -w 100 | sed 's/^/      /'
    eval ./configure "$PETSC_OPTS" \
        || { echo "ERROR: petsc configure failed; see configure.log" >&2; exit 4; }

    echo; echo "=== building petsc (make -j$NPROC) ==="
    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" -j"$NPROC" all \
        || { echo "ERROR: petsc build failed" >&2; exit 4; }
else
    echo "    petsc already built at $PETSC_DIR/$PETSC_ARCH, skipping"
fi

# --- verify the partitioners BEFORE going further ------------------------
# This is the whole reason for the build. If it fails here, stop: a Firedrake
# on top of a partitioner-less PETSc is exactly what we already have.
echo; echo "=== checking petsc has the parallel partitioners ==="
CONF_H="$PETSC_DIR/$PETSC_ARCH/include/petscconf.h"
ok=1
for k in PARMETIS PTSCOTCH; do
    if grep -q "define PETSC_HAVE_${k} 1" "$CONF_H" 2>/dev/null; then
        echo "    PETSC_HAVE_${k}: yes"
    else
        echo "    PETSC_HAVE_${k}: NO" >&2; ok=0
    fi
done
[ "$ok" = 1 ] || { echo "ERROR: partitioners missing; not proceeding" >&2; exit 5; }

# --- Firedrake -----------------------------------------------------------
cd "$PREFIX"
if [ ! -d "$VENV" ]; then
    echo; echo "=== creating venv ==="
    python3 -m venv "$VENV" || exit 6
fi
# shellcheck disable=SC1091
. "$VENV/bin/activate"

echo; echo "=== pip install firedrake ==="
python3 -m pip install --quiet --upgrade pip setuptools wheel
# firedrake-configure --show-env emits PETSC_DIR/PETSC_ARCH/HDF5 hints.
set -a; eval "$($FDC --show-env)"; set +a
python3 -m pip install --no-binary h5py 'firedrake[check]' \
    || { echo "ERROR: pip install firedrake failed" >&2; exit 6; }

# --- prove it works, in parallel, with a real partition ------------------
echo; echo "=== verifying ==="
export OMP_NUM_THREADS=1
export PYOP2_CACHE_DIR="$PREFIX/.pyop2_cache_build"
mkdir -p "$PYOP2_CACHE_DIR"

python3 -c "import firedrake; print('    firedrake', firedrake.__version__)" 2>&1 | tail -1

cat > "$PREFIX/_verify_part.py" <<'PYEOF'
from firedrake import *
from firedrake.petsc import PETSc
import sys
PETSc.Options()["petscpartitioner_type"] = sys.argv[1]
m = UnitSquareMesh(80, 80, distribution_parameters={"partition": True})
V = FunctionSpace(m, "CG", 1)
owned, total = V.dof_dset.size, V.dof_dset.total_size
g = m.comm.allreduce(total - owned) / max(m.comm.allreduce(owned), 1)
if m.comm.rank == 0:
    print(f"    {sys.argv[1]:9s} ghost/owned = {g:.4f}")
PYEOF

for p in parmetis ptscotch simple; do
    srun -n 4 python3 "$PREFIX/_verify_part.py" "$p" 2>/dev/null | tail -1 \
        || echo "    $p FAILED"
done

echo
echo "=== $(date -Is)  build done ==="
echo "Set these in site_env.sh:"
echo "    ISMIP7_FIREDRAKE=$VENV/bin/activate"
echo "    ISMIP7_MODULES=\"foss/2023b Python/3.11.5 CMake/3.27.6 M4/1.4.19 flex/2.6.4 Bison/3.8.2 libevent/2.1.12 zstd/1.5.5 Szip/2.1.1 libaec/1.0.6\""
