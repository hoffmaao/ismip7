# Slurm batch runners

Scheduler job scripts for the ISMIP7 icepack2 runs, on any cluster. Every wall
time and memory figure here was measured.

## How it works

| | |
|---|---|
| `sites/<cluster>.sh` | everything that differs between clusters: Firedrake venv, module loads, partitions by job length, account, node constraint, and paths that hold a checkout and 300 GB of forcing. |
| `site_env.sh` | picks the site file (`ISMIP7_SITE`, else hostname), refuses to run while a required field is empty, and provides `ismip7_activate`, `ismip7_banner` and the shared model defaults. |
| `submit.sh` | composes the scheduler command from the site file, so `submit.sh inversion ISMIP7_LC=2000` is the same line everywhere. |

The job scripts carry log paths and nothing else. A `#SBATCH` line is parsed
before any shell runs, so it cannot read a site file, and a header that
disagrees with the command line is rejected outright ("Requested node
configuration is not available"). A job submitted by hand without resources
stops immediately and says so.

```bash
submit.sh inversion  ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_5000_2000_buffered0.msh
submit.sh projection ISMIP7_EXPERIMENT=ssp585_cesm_waccm ISMIP7_OUTPUT=1
submit.sh smoke                                  # minutes, debug partition
submit.sh inversion --dry-run                    # print the sbatch line only
submit.sh projection --partition debug --time 00:30:00 --tasks 8
```

Any `KEY=VALUE` argument is exported into the job, which is how a run is
configured (knobs in `antarctica/README.md`). `--tasks`, `--mem`, `--time`,
`--partition`, `--constraint`, `--account` and `--name` override site defaults
for one submission, in either `--opt value` or `--opt=value` form.

## Adding your cluster

```bash
cp sites/template.sh sites/mycluster.sh
$EDITOR sites/mycluster.sh          # required fields are marked
ISMIP7_SITE=mycluster submit.sh smoke --dry-run
```

`site_recon.sh` on a login node prints what the file needs: associations and
accounts, partitions and their limits, node features and memory, queue depth.
It only reads. Add your hostname pattern to `ISMIP7_SITE_MATCH` and
`ISMIP7_SITE` becomes unnecessary.

Required: `ISMIP7_FIREDRAKE`, `ISMIP7_PART_LONG`, `ISMIP7_PART_SHORT`,
`ISMIP7_PART_DEBUG`, `ISMIP7_REPO`, `ISMIP7_WORK`. Everything else defaults.
A missing value is reported at submission with the file and variable named.

Job sizes come in pairs so inversions and forwards can differ:
`ISMIP7_TASKS_INV` and `ISMIP7_MEM_INV`, `ISMIP7_TASKS_FWD` and
`ISMIP7_MEM_FWD`, each falling back to `ISMIP7_TASKS` and `ISMIP7_MEM`. The
node constraint splits the same way (`ISMIP7_CONSTRAINT_INV`, `_FWD`).

### The sites that ship

| file | state |
|---|---|
| `sites/rice_nots.sh` | complete, in production. Details below. |
| `sites/iu_quartz.sh` | complete, from David Lilien's Quartz runners on the upstream `timing_matrix` branch: partition `general` (`debug` for tests), account `r00905`, the IU module stack (`module use /N/u/dlilien/Quartz/modulefiles`, then gnu, openmpi, python, zlib, hdf5, openblas, patchelf, petsc, firedrake), 16 ranks per node. A second IU user points `ISMIP7_FIREDRAKE` at their own build and `ISMIP7_ACCOUNT` at their own allocation. |
| `sites/uchicago_midway.sh` | a stub. Nobody has run this pipeline at RCC, so the required fields are empty and the first submission refuses until they are filled from RCC's documentation and `sinfo -s`. |
| `sites/local.sh` | no scheduler: every setting comes from the environment. For debugging a job script on a workstation, and what the chain tests use. |

## Rice NOTS, in detail

Measured on the login node in September 2026. `sites/rice_nots.sh` is the file
it produced.

The default account is `commons`, so `-A` is unnecessary. Four partitions are
granted:

| partition | wall limit | cap | nodes | cascadelake nodes |
|---|---|---|---|---|
| commons | 1 day | 1536 cpu / 2304 GB | 142 | 72 |
| long | 3 days | 384 cpu / 1 TB | 61 | 1 |
| debug | 30 min | 1 job | 2 | |
| scavenge | 1 hour | preemptible | 154 | many |

Cascade Lake is reachable today on `commons` and `scavenge`. Pinning it on
`long` would queue behind a single machine, so `sites/rice_nots.sh` sets
`ISMIP7_CONSTRAINT_INV=sapphirerapids` (`long` is where the 2 km inversions
run) and `ISMIP7_CONSTRAINT_FWD=cascadelake`, the generation the timings below
were measured on and the one the Firedrake build matches.

**deepsC is refused.** `sbatch --test-only -p deepsC` returns "Invalid account
or account/partition combination". It needs the separate `deepsc` account,
described in Slurm as Prof. Alan Levander's, with 76 users already in it and no
Slurm coordinators, so adding a user goes through the help desk. It is worth
having because it has no wall-time limit (a 2500 m inversion runs 20 to 26
hours against the one-day `commons` limit) and is less contended (182 pending
against 624). It does not solve memory: all 52 deepsC nodes are 187 GB, so the
2 km / 5 km-interior inversion at about 255 GB needs `long` and its 515 GB
Sapphire Rapids nodes.

Once granted, switch on the wrapper. Clear the constraint as well, since every
deepsC node is Cascade Lake while inversions default to Sapphire Rapids:

```bash
submit.sh inversion --partition deepsC --constraint '' --time 2-00:00:00
```

deepsC nodes are uniform: 40 physical cores (2 x 20), 80 logical
(`ThreadsPerCore=2`), 187135 MB, features `cascadelake,opath`, no wall limit.
The big-memory nodes live elsewhere (`commons` and `long` reach 1.5 TB, `as143`
748 GB). The solver is memory-bandwidth bound, so `submit.sh` passes
`--hint=nomultithread` and ranks land on physical cores. The CPU matches the
workstation the timings came from (Xeon Gold 5218R, 2.10 GHz), so they transfer
one to one and are conservative, the workstation having been oversubscribed.

## Building Firedrake on a cluster

The lessons generalise. The module names do not.
`build_firedrake_rice.sbatch` is the worked example and runs at Rice only.
Another site copies it to `build_firedrake_<site>.sbatch` and substitutes its
own stack; `submit.sh build` says so when the site is not `rice_nots`.

Five submissions were needed at Rice, each failure a gap between login and
compute nodes:

| what broke | why | fix |
|---|---|---|
| `firedrake-configure` refused | RHEL 9.4 is unsupported | `--os unknown`, PETSc builds deps from source |
| `git: command not found` | compute nodes have no git | fetch the PETSc tarball with `curl` |
| Bison would not build | no `flex`, `bison`, `m4` | load `M4/1.4.19 flex/2.6.4 Bison/3.8.2` |
| PnetCDF configure failed | OpenMPI adds `-levent_core`, libevent off `LIBRARY_PATH` | load `libevent/2.1.12` |
| NetCDF link failed | `ld: cannot find -lzstd` | load `zstd/1.5.5 Szip/2.1.1 libaec/1.0.6` |
| modules "loaded" but gcc was 11 | a bare name like `flex` resolves against the wrong GCCcore and drops the set back to `/usr/bin` | pin every version, abort if `gcc` is not under `/opt/apps` |

The last one is silent: `module load` returns 0 and you get gcc 11, Python 3.9
and no `mpicc`. The script fails loudly now. While debugging, never pipe
`module load` through `head` or `grep`; SIGPIPE breaks the load and you measure
an environment the job never sees.

On partitioners: upstream's own option set carries `--download-ptscotch`, so a
stock build works. The workstation build was configured `--with-ptscotch`
against a system scotch that supplied nothing, and PETSc fell back to `simple`,
which cuts by point index. That is harmless on a structured mesh and
catastrophic on a gmsh-ordered one (ghost/owned 16.0 at 4 ranks, 287 at 32).
The build adds `--download-parmetis` so two partitioners can be cross-checked.

Rice has no Firedrake module, so the build uses EasyBuild `foss/2023b` (GCC
13.2.0, OpenMPI 4.1.6, OpenBLAS, ScaLAPACK, FFTW) with Python 3.11.5 on the
same GCCcore, CMake 3.27.6, and the M4, flex, Bison, libevent, zstd, Szip,
libaec and libGLU modules the gaps above require. `ISMIP7_MODULES` must name
the exact set the venv was built against, which a successful build prints at
the end.

### After the build

`install_deps.sh` (once, on a login node) pip-installs the data stack and the
four editable packages: `icepack`, `icepack2`, `tlm_adjoint`, `icepack_tools`.
It takes the venv and work filesystem from the site file and expects the
sources under `$ISMIP7_WORK/sw/src`. They are rsynced from a workstation rather
than cloned, since `icepack2` carries uncommitted edits the inversion needs and
`icepack_tools` is private. Two more gaps surfaced here: `/tmp` is not writable
on the login nodes (the script sets `TMPDIR`), and the `gmsh` wheel dlopens
`libGLU.so.1`.

`verify.sbatch` proves the build works across ranks: four tasks under `srun`,
each partitioner on a unit square, then the real 2500 m mesh under `ptscotch`
and `simple`. The build job's own check ran `srun -n 4` inside a one-task
allocation, which fails for every type and proves nothing.

`smoke.sbatch` runs the 32 km inversion for two iterates on four ranks, a few
minutes on `scavenge`, so a multi-day job cannot die on a missing file hours
after queueing. Run it after any change to the stack.

## The job scripts

### `partition_probe.sbatch`, run first after a build

Distributes the mesh at 1 to 32 ranks and reports the ghost-to-owned dof ratio.
Hundredths mean a working parallel partitioner. Tens or hundreds mean PETSc
fell back to `simple` and every scaling number would measure that. The
workstation build reports 2.97 at 2 ranks and 287 at 32, which is why no
scaling curve is quoted from it.

### `inversion.sbatch`, self-resuming

Defaults write `inversion_icepack2_rc_n3_dg0_logvelnet_2500.h5` under the
settings the 2500 m result came from: sigma-normalised velocity misfit with
ISSM's logarithmic term, the pointwise dH/dt term, and the integrated net
mass-balance constraint that is off by default in the repo.

`site_env.sh` defaults `ISMIP7_FRICTION` to `regularized_coulomb` everywhere.
Budd's shelf gate was a sign test on the roundoff residue of the effective
pressure, so Budd MAPs predating the fix need re-inverting; pass
`ISMIP7_FRICTION=budd` for those. The law picks the filename tag, so a Budd
re-inversion writes its own `_budd` MAP.

`submit_inversions.sh [B|C|BC]` submits the two 2 km strategies (B: cell-mean
BedMachine sampling on the 20 km-interior mesh; C: the 5 km-interior mesh with
vertex sampling), each named so the converged 2500 m MAP is never touched, with
partition, constraint and rank count from the site file.

**The chain.** A 2 km inversion can outlast a wall limit. The inversion
checkpoints `ISMIP7_MAP_OUT` every 20 iterates and warm-starts from it, so each
job queues its successor first with `--dependency=afterany`, copying partition,
constraint, memory, time and task layout from `scontrol`. Every link exits at
once if `<map>.done` exists, meaning the MAP reached disk. The driver writes
that marker as soon as the checkpoint write returns, so a kill in the tail
(final solve, summary figure) cannot lose it; the runner's post-`srun` grep for
the driver's `Saved MAP:` line is the fallback. Depth is capped by
`ISMIP7_CHAIN_MAX` (4). A warm start whose mesh dof ordering differs from the
run's own is refused, since a rank-count change mid-chain would scramble theta
and phi silently. Regression test: `tests/test_inversion_chain.py`.

### `projection.sbatch`, self-chaining

A 285-year projection is about five days at 2500 m, so this resubmits itself
with `--dependency=afterok` until the run reaches its end year, resuming
through `ISMIP7_AUTO_RESUME=1`. Each driver owns its end year and the runner
does not default `ISMIP7_T_END`; the chain reads the value the run used from
the driver's `Time-stepping: <start>-><end>` line. 24 h buys roughly 55
simulated years, so a full projection is about six links.

```bash
submit.sh projection ISMIP7_EXPERIMENT=control
```

| `ISMIP7_EXPERIMENT` | core | period |
|-------|------|----------------|
| `control` | 9 or 10 (by `ISMIP7_ESM`) | 2015-2300 |
| `ssp126_cesm_waccm` / `ssp126_mri_esm2` | 5 / 6 | 2015-2300 |
| `ssp370_cesm_waccm` / `ssp370_mri_esm2` | 3 / 4 | 2015-2100 |
| `ssp585_cesm_waccm` / `ssp585_mri_esm2` | 7 / 8 | 2015-2300 |
| `ocx` | 11 | 1979-2025 |
| `hist_cesm_waccm` / `hist_mri_esm2` | 1 / 2 | 1850-2014 |

**The chain stops on a non-zero exit and never retries.** The July
grounding-line blow-up looked like a run that needed more time, and chaining
through it would have burned days.

Each job asks its own final checkpoint what happened, and that verdict owns the
exit code. `stalled=1` (rescue ladder and subcycles both exhausted) exits 1 and
submits nothing, so a stall is never mailed out as a completed job. A clean
exit at the end year finishes. A missing `Time-stepping:` line stops the chain,
since the job cannot tell where it got to. A clean exit short of the end year
resubmits, unless `ISMIP7_CHAIN=0`, or auto-resume is off (a successor would
cold-start and repeat the years), or the final checkpoint is unreadable, or the
job advanced no years at all (setup alone spent the budget).

A relaunch from a stalled state has cleared the diagnostic-Newton wall in every
observed case (see Known issues in `antarctica/README.md`). The chain leaves
that judgement to you: resubmit, and auto-resume picks the run up.

The wall budget is derived per job. Each link reads its own partition's
`TimeLimit`, holds back 25 minutes and passes the rest as
`ISMIP7_WALL_STOP_MIN`, so the model stops a step early and writes a complete
final checkpoint for the successor.

## Measured costs

| run | ranks | memory | wall | core-hours |
|---|---|---|---|---|
| 2500 m inversion, 60-80 iterates | 12 | 80 GB | 1 day | 300 |
| 2500 m forward, per simulated year | 12 | 80 GB | 26 min | 5 |
| 2500 m projection, 285 years | 12 | 80 GB | 5 days | 1,470 |
| 2 km / 20 km interior inversion | 12 | 120 GB | 1.5 days | 450 |
| 2 km / 5 km interior inversion | 12 | 255 GB | 2-3 days | 800 |
| full 11-experiment set at 2500 m | | | | 15,000 |

Per iterate at 2500 m on 12 ranks: forward median 1081 s (p10 932, p90 1365),
adjoint 92 s, iterate 1174 s. The adjoint is 8% of the iterate, so the cost
sits in the forward's Newton continuation. Per forward step: 155 s median, of
which the level set is 13 s.

The eleven-experiment set is about nineteen node-days, and the experiments are
independent, so a handful of nodes finishes it inside a week.

## Open items

- Rank count is fixed at 12 because that is what was measured. Going higher is
  meaningful once the partition probe comes back clean.
- Every solve is a direct MUMPS factorisation, which sets the memory. The fix
  if large-memory nodes get scarce is a field split that eliminates the
  cell-wise stress and traction blocks and puts multigrid on the velocity
  operator.
