r"""staging.py - read a checkpoint from node-local disk instead of over NFS.

Every rank reading a checkpoint over NFS through ROMIO (``OMPI_MCA_io=romio321``,
chosen in ``site_core.sh`` for its write speed) sits in NFS lock round trips:
on NOTS /scratch, 30 Sep 2026, 32 ranks took 45+ min to load a 580 MB 2 km
MAP, about 12 KB per rank per 20 s. A local copy takes a second, so a run
that loads a MAP or a restart checkpoint at setup copies it first.

The copy is made once per node by that node's first rank and opened by every
rank. Only a single-node job stages: MPI-IO opens one file on all ranks, and
copies on two nodes' disks are two files.
"""
import atexit
import os
import shutil
import tempfile

from mpi4py import MPI

_COPIES = {}


def staging_enabled():
    r"""``ISMIP7_STAGE_READS``: ``1`` stages, ``0`` reads in place. The
    default stages under Slurm (a networked results tree) and reads in place
    off it (a workstation's own disk)."""
    value = os.environ.get("ISMIP7_STAGE_READS")
    if value is None:
        return bool(os.environ.get("SLURM_JOB_ID"))
    return value.strip().lower() not in ("", "0", "false", "no", "off")


def node_local_copy(path, comm=MPI.COMM_WORLD):
    r"""``path``, or a copy of it on this node's local disk.

    Collective over ``comm``. Returns ``path`` itself when staging is off,
    when ``comm`` spans more than one node, or when the copy fails on any
    node; otherwise every rank gets the same local path, under ``$TMPDIR``
    (``/tmp`` when unset), removed when the process exits. A path staged
    once is not copied again.
    """
    if path in _COPIES:
        return _COPIES[path]
    if not staging_enabled():
        return path
    node = comm.Split_type(MPI.COMM_TYPE_SHARED)
    try:
        if node.size != comm.size:
            return path
        local, ok = None, True
        if node.rank == 0:
            try:
                base = os.environ.get("TMPDIR") or "/tmp"
                tmpdir = tempfile.mkdtemp(prefix="ismip7_stage_", dir=base)
                atexit.register(shutil.rmtree, tmpdir, True)
                local = os.path.join(tmpdir, os.path.basename(path))
                shutil.copyfile(path, local)
            except OSError:
                ok = False
        local, ok = node.bcast((local, ok), root=0)
        if not comm.allreduce(ok, op=MPI.LAND):
            return path
        _COPIES[path] = local
        return local
    finally:
        node.Free()
