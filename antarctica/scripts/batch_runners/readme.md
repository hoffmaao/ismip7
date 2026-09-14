# Slurm batch runners

Scheduler job scripts for the ISMIP7 icepack2 runs, on any cluster. Every wall
time and memory figure below was measured, not estimated.

## How it works

Three pieces, and the split matters:

| | |
|---|---|
| `sites/<cluster>.sh` | everything that differs between clusters: the Firedrake venv, the module loads it was built against, partitions by job length, the account to charge, a node constraint, and the paths that can hold a checkout and 300 GB of forcing. |
| `site_env.sh` | picks the site file (from `ISMIP7_SITE`, else by hostname), refuses to run while a required field is empty, and turns it into an environment: `ismip7_activate`, `ismip7_banner`, and the model defaults every run shares. |
| `submit.sh` | composes the scheduler command from the site file. `submit.sh inversion ISMIP7_LC=2000` is the same line at Rice and at IU. |

The job scripts themselves (`inversion.sbatch`, `projection.sbatch`,
`smoke.sbatch`, `verify.sbatch`, `partition_probe.sbatch`,
`build_firedrake.sbatch`) carry **no resource directives at all**, only their
log paths. That is deliberate: a `#SBATCH` line is parsed before any shell
runs, so it cannot read a site file, and a header that disagrees with the
command line is not overridden but rejected outright ("Requested node
configuration is not available", which cost two submissions before the split).
A job submitted by hand without resources stops immediately and says so.

```bash
# the usual calls
submit.sh inversion  ISMIP7_LC=2000 ISMIP7_LC_COARSE=5000 ISMIP7_MESH=$PWD/antarctica/mesh/antarctica_5000_2000_buffered0.msh
submit.sh projection ISMIP7_EXPERIMENT=ssp585_cesm_waccm ISMIP7_OUTPUT=1
submit.sh smoke                                   # a few minutes, debug partition
submit.sh inversion --dry-run                     # print the sbatch line, submit nothing
submit.sh projection --partition debug --time 00:30:00 --tasks 8   # override anything
```

Any `KEY=VALUE` argument is exported into the job; that is how a run is
configured (the knobs are in `antarctica/README.md`). `--tasks`, `--mem`,
`--time`, `--partition`, `--constraint`, `--account` and `--name` override the
site defaults for one submission.

Inversions and forwards are sized separately: a site file may set
`ISMIP7_TASKS_INV`/`ISMIP7_MEM_INV` and `ISMIP7_TASKS_FWD`/`ISMIP7_MEM_FWD`,
each falling back to the single `ISMIP7_TASKS`/`ISMIP7_MEM`. At Rice the
forward runs at the 12 ranks and 96 GB it was measured at while the inversion
keeps 32 ranks and 240 GB.

## Adding your cluster

```bash
cp sites/template.sh sites/mycluster.sh
$EDITOR sites/mycluster.sh          # the REQUIRED fields are marked
ISMIP7_SITE=mycluster submit.sh smoke --dry-run
```

`site_recon.sh`, run on a login node, prints most of what the file needs: the
Slurm associations and partition limits, the node hardware, the queue depth.
It only reads, and submits nothing. Add your hostname pattern to the file's
`ISMIP7_SITE_MATCH` and `ISMIP7_SITE` stops being necessary.

Required, because a job cannot start without them: `ISMIP7_FIREDRAKE`,
`ISMIP7_PART_LONG`, `ISMIP7_PART_SHORT`, `ISMIP7_PART_DEBUG`, `ISMIP7_REPO`,
`ISMIP7_WORK`.
Everything else has a working default. A missing value is reported at
submission with the file and the variable named, never minutes into a queued
job.

### The sites that ship

| file | state |
|---|---|
| `sites/rice_nots.sh` | complete and in production (the measured details below). |
| `sites/iu_quartz.sh` | complete, taken from David Lilien's own Quartz runners on the upstream `timing_matrix` branch: partition `general` (`debug` for tests), account `r00905`, the IU module stack (`module use /N/u/dlilien/Quartz/modulefiles`, then gnu/openmpi/python/zlib/hdf5/openblas/patchelf/petsc/firedrake), 16 ranks per node. A second IU user changes `ISMIP7_FIREDRAKE` to their own build and `ISMIP7_ACCOUNT` to their own allocation. |
| `sites/uchicago_midway.sh` | **a stub**: nobody has run this pipeline at RCC yet, so the required fields are empty rather than guessed and the first submission will refuse until they are filled from RCC's documentation and `sinfo -s`. |
| `sites/local.sh` | no scheduler at all: takes every setting from the environment. For debugging a job script on a workstation, and what the chain tests use. |

## Rice NOTS, in detail

The reference site: everything in this section was measured on the login node
in September 2026, and `sites/rice_nots.sh` is the file it produced.

### Slurm access at Rice: resolved, except deepsC

RESOLVED 2026-09-12. The association exists and the default account is
`commons`, so `-A` is unnecessary. A `--test-only` submit succeeds on four
partitions:

| partition | wall limit | your cap | nodes |
|---|---|---|---|
| commons | 1 day | 1536 cpu / 2304 GB | 142 |
| long | 3 days | 384 cpu / 1 TB | 61 |
| debug | 30 min | 1 job | 2 |
| scavenge | 1 hour | preemptible | 154 |

**You already have Cascade Lake.** `commons` carries 72 Cascade Lake nodes of
142, and scavenge has plenty, so `--constraint cascadelake` reaches
EEPS-generation hardware today. Do NOT pin that constraint on `long`: it has
exactly ONE Cascade Lake node, so the job would queue behind a single machine.
That is why `sites/rice_nots.sh` pins `sapphirerapids` instead, since `long` is
the partition the 2 km inversions need.

| partition | cascadelake nodes | wall limit |
|---|---|---|
| deepsC | 52 of 52 | infinite |
| commons | 72 of 142 | 1 day |
| long | **1** of 61 | 3 days |

**deepsC is still refused**, and only deepsC:

```
$ sbatch --test-only -p deepsC -N1 -n1 -t 00:01:00 --wrap=true
allocation failure: Invalid account or account/partition combination specified
```

It needs the separate `deepsc` account. The Slurm description names it
"alan lavendar`s deepsc account", i.e. Prof. Alan Levander, and 76 users are
already in it. No Slurm coordinators are set, so adding a user goes through the
help desk, presumably with the owner's approval.

**Why it is worth having**, given Cascade Lake is already reachable on commons:
deepsC has NO wall-time limit. A 2500 m inversion runs 20 to 26 hours, which is
uncomfortably close to the one-day limit on commons, and the projections are far
longer. deepsC is also less contended (182 pending vs 624 on commons).

**What it does NOT solve:** memory. Every deepsC node is 187 GB, so the 2 km /
5 km-interior inversion at ~255 GB still cannot run there. That one needs
`long`, where the 515 GB Sapphire Rapids nodes live. Those are newer and wider
than Cascade Lake (192 cpus), so they are better hardware for getting that run
done; the only cost is that timings taken there are not comparable with the
Cascade Lake numbers in this file.

Switching once granted is three options on the wrapper, which override the site
file for that submission. The constraint has to be cleared as well as the
partition: `sites/rice_nots.sh` pins `sapphirerapids`, every deepsC node is
Cascade Lake, and the two together can never be satisfied.

```
antarctica/scripts/batch_runners/submit.sh inversion \
    --partition deepsC --constraint '' --time 2-00:00:00
```

### What deepsC actually is (measured, and it corrects an earlier claim)

Every deepsC node is identical:

| | value |
|---|---|
| nodes | 52 |
| physical cores | 40 (2 sockets x 20) |
| logical CPUs | 80 (`ThreadsPerCore=2`) |
| memory | 187135 MB, about 182 GiB usable |
| features | `cascadelake,opath` |
| wall-time limit | **infinite** |

**An earlier version of these notes claimed deepsC had 768 GB and 1.5 TB tiers.
That was wrong.** Those nodes exist on NOTS, but in other partitions:

| partition | nodes | memory | limit |
|---|---|---|---|
| deepsC | 52 | 187 GB, uniform | infinite |
| commons | 142 | 187 GB and up, incl. 1.5 TB | 1 day |
| long | 61 | 257 GB and up, incl. 1.5 TB | 3 days |
| as143 | 8 | 748 GB | infinite |
| scavenge | 154 | 187 GB and up | 1 hour |

The consequence is concrete: **the 2 km / 5 km-interior inversion at ~255 GB
cannot run on deepsC.** It needs `long`, which caps at 3 days. The 20 km-interior
2 km inversion at ~120 GB does fit deepsC.

The cores are hyperthreaded. 40 are physical; Slurm advertises 80. The solver is
memory-bandwidth bound, so `submit.sh` passes `--hint=nomultithread` and places
ranks on physical cores. Treating 80 as usable would halve per-rank bandwidth.

The CPU is the same generation and clock as the workstation the timings were
measured on (Xeon Gold 5218R, Cascade Lake, 2.10 GHz), so they transfer nearly
one to one, and they are conservative: the workstation was oversubscribed
throughout, so a dedicated node should beat them.

## Building Firedrake on a cluster

The lessons generalise; the module names do not. Worked example: Rice NOTS.

### What compute nodes tend NOT to have

`build_firedrake.sbatch` builds PETSc + Firedrake once into
`/projects/ah301/sw`. Five submissions were needed to get it running, and every
failure was a gap between the login node and the compute nodes, or between a
module name and what it actually resolves to. They are recorded here because
each one costs a queue round trip to rediscover.

| what broke | why | fix in the script |
|---|---|---|
| `firedrake-configure` refused to run | NOTS is RHEL 9.4, not a supported OS | `--os unknown`, so PETSc builds deps from source |
| `git: command not found` | compute nodes have **no git at all**; the login node does | fetch the PETSc release tarball with `curl` |
| Bison would not build | no `flex`, `bison` or `m4` on a compute node | load `M4/1.4.19 flex/2.6.4 Bison/3.8.2` |
| PnetCDF configure failed | OpenMPI's wrapper adds `-levent_core` but libevent is not on `LIBRARY_PATH` | load `libevent/2.1.12` |
| NetCDF failed to link | same again: `ld: cannot find -lzstd` | load `zstd/1.5.5 Szip/2.1.1 libaec/1.0.6` |
| modules silently "loaded" but gcc was 11 | a bare name like `flex` resolves against the wrong GCCcore and leaves M4 unloaded, dropping the whole set back to `/usr/bin` | pin every version, then abort if `gcc` is not under `/opt/apps` |

The last one is the nastiest, because nothing errors. `module load` returns 0 and
you get system gcc 11, Python 3.9 and no `mpicc`. The script now fails loudly
instead. A related trap while debugging: piping `module load` through `head` or
`grep` breaks the load via SIGPIPE, so the environment you then measure is not
the one the job gets. Redirect to a file instead.

The partitioner question that motivated all of this turned out to have a
simpler answer than expected. Upstream's own option set already carries
`--download-ptscotch`, so a stock Firedrake build fixes it. The workstation
build is broken because it was configured `--with-ptscotch` against a system
scotch that never supplied the parallel partitioner, and PETSc silently fell
back to `simple`. Forcing each type through a real 4-rank distribution on the
workstation:

```
parmetis  -> PETSc error 56 (backend absent)
ptscotch  -> PETSc error 56
chaco     -> PETSc error 56
simple    -> works, ghost/owned 0.11 on a structured mesh
```

`simple` cuts by point index, which is harmless on a structured mesh and
catastrophic on a gmsh-ordered one: the Antarctica mesh measures ghost/owned
16.0 at 4 ranks rising to 287 at 32. The build adds `--download-parmetis` on
top so there are two partitioners to cross-check.

### After the build: the rest of the stack, and a smoke test

`install_deps.sh` (run once on a login node) pip-installs the data stack
and the four editable packages the inversion imports: `icepack`, `icepack2`,
`tlm_adjoint` and `icepack_tools`. They are rsynced from the workstation into
`/projects/ah301/sw/src`, not cloned, because `icepack2` carries two
uncommitted edits the inversion depends on and `icepack_tools` has no remote.
Two more compute-node gaps surfaced here: `/tmp` is not writable on the login
nodes (the script sets `TMPDIR`), and the `gmsh` wheel dlopens `libGLU.so.1`,
which the nodes lack, so `libGLU/9.0.3` joined `ISMIP7_MODULES`. HDF5 versions
agree (h5py 1.14.6 built against PETSc's, netCDF4 wheel 1.14.6).

`verify.sbatch` then proves the build works ACROSS RANKS: four tasks under
`srun`, distributing a unit square with each partitioner in turn and then the
real 2500 m mesh under `ptscotch` and `simple`. It exists because the build
job's own check ran `srun -n 4` inside a one-task allocation, which fails for
every type, `simple` included, and so says nothing about the build.

`smoke.sbatch` then runs the 32 km inversion for two iterates on four
ranks, a few minutes on `scavenge`, so that a multi-day job cannot die on a
missing file hours after it queued. Run it after any change to the stack.

### Toolchain for the Firedrake build

There is no Firedrake module on NOTS, so `build_firedrake.sbatch` builds
it against EasyBuild modules. The set it pins is `foss/2023b`, the toolchain
that actually resolves here (GCC 13.2.0, OpenMPI 4.1.6, OpenBLAS, ScaLAPACK,
FFTW) with Python 3.11.5 built against the same GCCcore so the ABI matches,
plus CMake 3.27.6 and the M4 / flex / Bison / libevent / zstd / Szip / libaec
set the compute-node gaps above make necessary. An earlier version of these
notes named `foss/2025b`; the build does not use it. `ISMIP7_MODULES` must name
the *exact* set the venv was built against, which the successful build prints
at the end, or MPI mismatches at runtime rather than at submission.

## The job scripts

### `partition_probe.sbatch` - run first, after a new build

Distributes the mesh at 1 to 32 ranks and reports the ghost-to-owned dof ratio.
No solve, no data, half an hour of queue. It answers the only question that
matters before any timing: does this build have a working parallel graph
partitioner? Hundredths mean yes. Tens or hundreds mean PETSc fell back to its
`simple` partitioner and every scaling number would measure that instead of the
model. The workstation build reports 2.97 at 2 ranks rising to 287 at 32, which
is 284x to 6330x the ideal, and is why no scaling curve has been quoted from it.

### `inversion.sbatch` - self-resuming

Defaults write `inversion_icepack2_rc_n3_dg0_logvelnet_2500.h5` under the
settings the 2500 m result came from: the sigma-normalised velocity misfit
with ISSM's logarithmic term, the pointwise dH/dt term, and the integrated
net mass-balance constraint that is off by default in the repo.

The friction law is the one deliberate difference from that run.
`site_env.sh` defaults `ISMIP7_FRICTION` to `regularized_coulomb` at every
site, because every inversion now runs that law. Budd's shelf gate was a sign test on the
roundoff residue of the effective pressure, so the Budd MAPs that predate the
fix have to be re-inverted; pass `ISMIP7_FRICTION=budd` explicitly for those.
The law also picks the filename tag, so a Budd re-inversion writes its own
`_budd` MAP instead of resuming from and then overwriting the RC one.

`submit_inversions.sh [B|C|BC]` submits the two 2 km
strategies (B: cell-mean BedMachine sampling on the 20 km-interior mesh; C: the
5 km-interior mesh with vertex sampling), each named so the converged 2500 m
map is never touched. Both take this cluster's inversion defaults through
`submit.sh`, so the partition, constraint and rank count come from the site
file. At Rice that is `long`, Sapphire Rapids and 32 ranks; at IU it is
`general`, no constraint and 16 ranks.

**The chain.** A 2 km inversion can outlast a wall limit, and the inversion
checkpoints `ISMIP7_MAP_OUT` every 20 iterates and can warm-start from it. So
each job queues its own successor FIRST with `--dependency=afterany`, copying
its partition, constraint, memory, time and task layout from `scontrol`; every
link exits immediately if `<map>.done` exists, which means the MAP reached
disk (`ISMIP7_MAXITER` is the run's budget, so reaching the iteration cap is
the normal end; the successor is for a job the wall clock killed before the
MAP was written). The inversion driver writes the marker itself the moment the
checkpoint write returns, so a kill inside the tail after it (final solve,
summary figure) cannot lose it; the runner's post-`srun` rule, which greps the
log for the driver's `Saved MAP:` line, is the fallback for a driver that could
not write it. Regression test: `tests/test_inversion_chain.py`. Depth is capped
by `ISMIP7_CHAIN_MAX` (4). The inversion itself now refuses a warm start whose
mesh dof ordering differs from its own (a rank-count change mid-chain would
otherwise scramble theta/phi silently).

It writes to an explicit `ISMIP7_MAP_OUT`. Without that, any run, including a
smoke test, writes the production filename and can silently replace a
converged map.

### `projection.sbatch` - self-chaining

A 285-year projection is about five days at 2500 m and fits no ordinary queue,
so this resubmits itself with `--dependency=afterok` until the run reaches its
end year, resuming through `ISMIP7_AUTO_RESUME=1` from its own newest
checkpoint. Each driver owns that end year and the runner does not default
`ISMIP7_T_END`; the chain reads the value the run actually used back out of
the `Time-stepping: <start>-><end>` line the driver prints. Setting
`ISMIP7_T_END` overrides every experiment, so leave it unset unless you mean
to. 24 h buys roughly 55 simulated years, so a full projection is about six
chained jobs.

```
antarctica/scripts/batch_runners/submit.sh projection ISMIP7_EXPERIMENT=control
```

`ISMIP7_EXPERIMENT` selects the driver, one of the ten cores:

| value | core | period covered |
|-------|------|----------------|
| `control` | 9 / 10 (by `ISMIP7_ESM`) | 2015-2300 |
| `ssp126_cesm_waccm` | 5 | 2015-2300 |
| `ssp126_mri_esm2` | 6 | 2015-2300 |
| `ssp370_cesm_waccm` | 3 | 2015-2100 |
| `ssp370_mri_esm2` | 4 | 2015-2100 |
| `ssp585_cesm_waccm` | 7 | 2015-2300 |
| `ssp585_mri_esm2` | 8 | 2015-2300 |
| `ocx` | 11 | 1979-2025 |
| `hist_cesm_waccm` | 1 | 1850-2014 |
| `hist_mri_esm2` | 2 | 1850-2014 |

**The chain stops on a non-zero exit and does not retry.** That is deliberate.
The July grounding-line blow-up looked exactly like a run that just needed more
time, and chaining through it would have burned days producing nothing.

Beyond that, each job asks its OWN final checkpoint what happened before it
decides anything, and that verdict owns the exit code. A checkpoint written
with `stalled=1` (the in-run rescue ladder and the subcycles were both
exhausted) exits 1 and submits nothing, whether or not chaining was enabled, so
a stall is never mailed out as a completed job. A clean exit that reached the
end year read back from the driver's `Time-stepping:` line finishes; if that
line is missing the job cannot tell, and stops without resubmitting. A clean
exit short of the end year resubmits, except that `ISMIP7_CHAIN=0` disables
resubmission entirely; auto-resume off
(`ISMIP7_AUTO_RESUME=0`) stops the chain, because a successor inherits no
`ISMIP7_RESTART` and would cold-start and repeat the same years; an unreadable
final checkpoint stops it, since nothing then says where to resume; and a job
that advanced no years at all stops with exit 1, since its setup alone spent
the whole wall budget and a successor under the same budget would not advance
either.

A stall stops the chain, not the science: a relaunch from the saved state has
cleared the diagnostic-Newton wall in every case observed so far (see "Known
issues" in `antarctica/README.md`). The chain will not make that call for you,
so resubmit this script by hand when you judge the wall worth another process;
auto-resume picks the run up from its saved state.

That wall budget is derived per job, not inherited: each link reads its own
partition's `TimeLimit`, holds back 25 minutes and passes the rest to the
driver as `ISMIP7_WALL_STOP_MIN` (documented in the env table of
`antarctica/README.md`). The model then stops a step early and writes its final
checkpoint, so the successor resumes from a complete state rather than from
whatever survived being killed at the limit.

## Measured costs

| run | ranks | memory | wall | core-hours |
|---|---|---|---|---|
| 2500 m inversion, 60-80 iterates | 12 | 80 GB | ~1 day | ~300 |
| 2500 m forward, per simulated year | 12 | 80 GB | ~26 min | ~5 |
| 2500 m projection, 285 years | 12 | 80 GB | ~5 days | ~1,470 |
| 2 km / 20 km interior inversion | 12 | ~120 GB | ~1.5 days | ~450 |
| 2 km / 5 km interior inversion | 12 | ~255 GB | ~2-3 days | ~800 |
| full 11-experiment set at 2500 m | | | | ~15,000 |

Per-iterate detail at 2500 m on 12 ranks: forward median 1081 s (p10 932, p90
1365), adjoint 92 s, iterate 1174 s. The adjoint is 8% of the iterate, so the
cost is the forward's Newton continuation, not the adjoint. Per forward step:
155 s median, of which the level set is 13 s and the rest is the warm Newton
solve.

The eleven-experiment set is about nineteen node-days. The experiments are
independent of one another, so on even a handful of the twenty condo nodes that
is under a week of wall time.

## Open items

- ~~Wall-time limit on `deepsC` is unknown~~ RESOLVED 2026-09-09: it is
  infinite. The 24 h in the scripts is now a self-imposed courtesy on a shared
  condo rather than a constraint, and the chaining matters on `commons` and
  `long` instead.
- The queue was nearly empty when measured (4 jobs, 2 running, 29 of 52 nodes
  only partially allocated), so the condo is not contended today.
- Rank count is fixed at 12 because that is what was measured. Going higher is
  only meaningful once the partition probe comes back clean.
- There is no iterative solver path in the model today: every solve is a direct
  MUMPS factorisation, which is what sets the memory. If the 768 GB nodes turn
  out to be scarce, the fix is a field split that eliminates the cell-wise
  stress and traction blocks and puts multigrid on the velocity operator.
