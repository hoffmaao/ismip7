# Slurm batch runners for Rice NOTS

Cascade Lake job scripts for the ISMIP7 icepack2 runs. Every wall time and
memory figure below was measured, not estimated; the provenance is in the
table at the end.

## Before the first submission

Three site-specific values must be filled into `nots_env.sh`. They are left
empty on purpose, because a wrong module name fails minutes into a queued job
rather than at submission.

| variable | what it is | status as of 2026-09-09 |
|---|---|---|
| `NOTS_ACCOUNT` | Slurm account to charge | **BLOCKED, needs a CRC admin.** See below. |
| `NOTS_FIREDRAKE` | the Firedrake venv `activate` on NOTS | `nots_build_firedrake.sbatch` builds it to `/projects/ah301/sw/venv-firedrake` |
| `NOTS_MODULES` | the module loads that venv was **built** against | the exact pinned set is printed at the end of a successful build |

`ISMIP7_REPO` now defaults to `/projects/ah301/ismip7`. That directory exists
and is writable. The forcing tree is 316 GB and the input data 8 GB, and
`/home` is a 10 TB NFS export shared by the whole cluster, so neither belongs
there. `/scratch` is 1.1 PB but 88% full and subject to purge.

### Slurm access: resolved, except deepsC

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
142, and scavenge has plenty, so the job scripts pin `-C cascadelake` and run on
EEPS-generation hardware today. Do NOT pin the constraint on `long`: it has
exactly ONE Cascade Lake node, so the job would queue behind a single machine.

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

Switching once granted is one flag, since a command-line option beats the
`#SBATCH` directive:

```
sbatch -p deepsC --time=2-00:00:00 antarctica/scripts/batch_runners/nots_inversion.sbatch
```

### Building Firedrake here: what the compute nodes do NOT have

`nots_build_firedrake.sbatch` builds PETSc + Firedrake once into
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

`nots_install_deps.sh` (run once on a login node) pip-installs the data stack
and the four editable packages the inversion imports: `icepack`, `icepack2`,
`tlm_adjoint` and `icepack_tools`. They are rsynced from the workstation into
`/projects/ah301/sw/src`, not cloned, because `icepack2` carries two
uncommitted edits the inversion depends on and `icepack_tools` has no remote.
Two more compute-node gaps surfaced here: `/tmp` is not writable on the login
nodes (the script sets `TMPDIR`), and the `gmsh` wheel dlopens `libGLU.so.1`,
which the nodes lack, so `libGLU/9.0.3` joined `NOTS_MODULES`. HDF5 versions
agree (h5py 1.14.6 built against PETSc's, netCDF4 wheel 1.14.6).

`nots_smoke.sbatch` then runs the 32 km inversion for two iterates on four
ranks, a few minutes on `scavenge`, so that a multi-day job cannot die on a
missing file hours after it queued. Run it after any change to the stack.

### Toolchain for the Firedrake build

Firedrake is not built on NOTS. EasyBuild modules are available and the newest
coherent toolchain is `foss/2025b`. Confirmed present: GCC through 15.1.0,
OpenMPI 4.1.6, Python 3.11.5, CMake 3.27.6. Whatever is chosen, `NOTS_MODULES`
must name the *exact* set the venv was built against, or MPI mismatches at
runtime rather than at submission.

## What deepsC actually is (measured, and it corrects an earlier claim)

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
memory-bandwidth bound, so the job scripts pass `--hint=nomultithread` and place
ranks on physical cores. Treating 80 as usable would halve per-rank bandwidth.

The CPU is the same generation and clock as the workstation the timings were
measured on (Xeon Gold 5218R, Cascade Lake, 2.10 GHz), so they transfer nearly
one to one, and they are conservative: the workstation was oversubscribed
throughout, so a dedicated node should beat them.

## The scripts

### 1. `nots_partition_probe.sbatch` — run first, after the build

Distributes the mesh at 1 to 32 ranks and reports the ghost-to-owned dof ratio.
No solve, no data, half an hour of queue. It answers the only question that
matters before any timing: does this build have a working parallel graph
partitioner? Hundredths mean yes. Tens or hundreds mean PETSc fell back to its
`simple` partitioner and every scaling number would measure that instead of the
model. The workstation build reports 2.97 at 2 ranks rising to 287 at 32, which
is 284x to 6330x the ideal, and is why no scaling curve has been quoted from it.

### 2. `nots_inversion.sbatch` — self-resuming

Defaults reproduce `inversion_icepack2_budd_n3_dg0_logvelnet_2500.h5`: the
sigma-normalised velocity misfit with ISSM's logarithmic term, the pointwise
dH/dt term, and the integrated net mass-balance constraint that is off by
default in the repo. `nots_submit_inversions.sh [B|C|BC]` submits the two 2 km
strategies (B: cell-mean BedMachine sampling on the 20 km-interior mesh; C: the
5 km-interior mesh with vertex sampling) to `long` on Sapphire Rapids at 32
ranks, each named so the converged 2500 m map is never touched.

**The chain.** A 2 km inversion can outlast a wall limit, and the inversion
checkpoints `ISMIP7_MAP_OUT` every 20 iterates and can warm-start from it. So
each job queues its own successor FIRST with `--dependency=afterany`, copying
its partition, constraint, memory, time and task layout from `scontrol`; every
link exits immediately if `<map>.done` exists, which is written only when
L-BFGS-B returns for a reason other than its iteration cap. Depth is capped by
`ISMIP7_CHAIN_MAX` (4). The inversion itself now refuses a warm start whose
mesh dof ordering differs from its own (a rank-count change mid-chain would
otherwise scramble theta/phi silently).

It writes to an explicit `ISMIP7_MAP_OUT`. Without that, any run, including a
smoke test, writes the production filename and can silently replace a
converged map.

### 3. `nots_projection.sbatch` — self-chaining

A 285-year projection is about five days at 2500 m and fits no ordinary queue,
so this resubmits itself with `--dependency=afterok` until the run reaches
`ISMIP7_T_END`, resuming through `ISMIP7_AUTO_RESUME=1` from its own newest
checkpoint. 24 h buys roughly 55 simulated years, so a full projection is about
six chained jobs.

```
sbatch --export=ALL,ISMIP7_EXPERIMENT=control \
       antarctica/scripts/batch_runners/nots_projection.sbatch
```

`ISMIP7_EXPERIMENT` selects the driver: `control`, `ssp126_cesm_waccm`, `ocx`,
`hist_cesm_waccm`, `hist_mri_esm2`.

**The chain stops on a non-zero exit and does not retry.** That is deliberate.
The July grounding-line blow-up looked exactly like a run that just needed more
time, and chaining through it would have burned days producing nothing.
`ISMIP7_CHAIN=0` disables resubmission entirely.

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
