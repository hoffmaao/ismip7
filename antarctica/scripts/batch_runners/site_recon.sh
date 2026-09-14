#!/bin/bash
# Run this ON a NOTS login node. It gathers the site-specific facts the
# site file still has blank (sites/template.sh lists them), plus partition limits.
# It only reads; it submits nothing and changes nothing.
echo "===== 1. slurm account (fills ISMIP7_ACCOUNT) ====="
sacctmgr -n show assoc user="$USER" format=account%30 2>&1 | sort -u

echo; echo "===== 2. deepsC limits (wall time, node count, state) ====="
sinfo -p deepsC -o "%20P %5a %12l %6D %6t %N" 2>&1

echo; echo "===== 3. deepsC hardware (cores, memory, feature names) ====="
sinfo -p deepsC -N -o "%20N %5c %10m %30f" 2>&1 | sort -u

echo; echo "===== 4. current queue depth ====="
squeue -p deepsC -o "%.10i %.20j %.10u %.2t %.11M %.6D %R" 2>&1 | head -25

echo; echo "===== 5. is Firedrake already built anywhere? (fills ISMIP7_FIREDRAKE) ====="
for d in "$HOME" /projects /scratch /work; do
    [ -d "$d" ] || continue
    find "$d" -maxdepth 5 -name activate -path '*firedrake*' 2>/dev/null | head -10
done
echo "(empty means Firedrake is not built yet on NOTS)"

echo; echo "===== 6. toolchain modules available (fills ISMIP7_MODULES) ====="
module -t avail 2>&1 | grep -iE '^(gcc|openmpi|mpich|intel|python|cmake|petsc)' | head -40

echo; echo "===== 7. storage quotas (the forcing tree is 316 GB) ====="
echo "-- home --";     du -sh "$HOME" 2>/dev/null
for v in SCRATCH WORK PROJECT; do
    eval "p=\$$v"; [ -n "$p" ] && echo "-- $v = $p --" && df -h "$p" 2>/dev/null | tail -1
done
quota -s 2>/dev/null || echo "(no quota command)"
