#!/bin/bash
# Run this on any Slurm login node. It prints what a site file needs, in the
# order sites/template.sh asks for it, so the output is what you paste into
# sites/<name>.sh. It only reads; it submits nothing and changes nothing.
#
#   bash antarctica/scripts/batch_runners/site_recon.sh | tee recon.txt

echo "===== 1. accounts you may charge (fills ISMIP7_ACCOUNT) ====="
sacctmgr -n show assoc user="$USER" format=account%30,partition%20 2>&1 | sort -u

echo; echo "===== 2. partitions and their limits (fills ISMIP7_PART_LONG/SHORT/DEBUG) ====="
echo "    pick the longest wall limit for LONG, a day or so for SHORT, the"
echo "    shortest for DEBUG."
sinfo -o "%20P %5a %12l %6D %10T" 2>&1 | sort -u

echo; echo "===== 3. node hardware and features (fills ISMIP7_CONSTRAINT, ISMIP7_TASKS, ISMIP7_MEM) ====="
echo "    %c is logical CPUs: halve it if the nodes are hyperthreaded, because"
echo "    the solver is memory-bandwidth bound and wants physical cores."
sinfo -N -o "%20P %5c %10m %30f" 2>&1 | sort -u

echo; echo "===== 4. queue depth, per partition ====="
squeue -h -o "%P" 2>&1 | sort | uniq -c | sort -rn | head -25

echo; echo "===== 5. is Firedrake already built anywhere? (fills ISMIP7_FIREDRAKE) ====="
for d in "$HOME" "${SCRATCH:-}" "${WORK:-}" "${PROJECT:-}" /projects /scratch /work; do
    [ -n "$d" ] && [ -d "$d" ] || continue
    find "$d" -maxdepth 5 -name activate -path '*firedrake*' 2>/dev/null | head -10
done
echo "(empty means Firedrake is not built here yet; see build_firedrake_rice.sbatch)"

echo; echo "===== 6. toolchain modules available (fills ISMIP7_MODULES) ====="
module -t avail 2>&1 | grep -iE '^(gcc|gnu|foss|openmpi|mpich|intel|python|cmake|petsc|hdf5|firedrake)' | head -40

echo; echo "===== 7. storage, for ISMIP7_REPO and ISMIP7_WORK ====="
echo "    the forcing tree is about 313 GB for everything, 25 GB per scenario."
echo "-- home --";     df -h "$HOME" 2>/dev/null | tail -1
for v in SCRATCH WORK PROJECT; do
    eval "p=\${$v:-}"; [ -n "$p" ] && echo "-- $v = $p --" && df -h "$p" 2>/dev/null | tail -1
done
quota -s 2>/dev/null || echo "(no quota command)"
