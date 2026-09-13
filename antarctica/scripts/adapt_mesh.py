r"""Adapt the mesh of a forward checkpoint the way Úa's AdaptMesh does.

    mpiexec -n 4 python antarctica/scripts/adapt_mesh.py CHECKPOINT.h5 \
        --out-checkpoint NEW.h5 [--out-mesh NEW.msh] [--rebuild-aref]

Reads the checkpoint's mesh and state, builds Úa's desired-element-size field
from the ISMIP7_ADAPT_* configuration (icepack2_tools.adapt_mesh), remeshes
the domain globally with gmsh on rank 0, transfers the state onto the new mesh
(MapFbetweenMeshes semantics: interpolation, ThickMin outside, bed re-sampled
from BedMachine, surface by flotation) and writes a checkpoint the forward can
restart from. The new mesh and its boundary-id sidecar land next to the old
mesh in antarctica/mesh/.

--rebuild-aref drops the frozen apparent-mass-balance reference so the forward
rebuilds it on the new mesh; only legitimate for a t=0 state (the initial
adaptation), never mid-run.
"""
import argparse
import glob
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

import firedrake as fd  # noqa: E402
import rasterio  # noqa: E402
from firedrake import Function, FunctionSpace  # noqa: E402
from firedrake.petsc import PETSc  # noqa: E402

from icepack2_tools.adapt_mesh import (AdaptMeshConfig, desired_element_size,  # noqa: E402
                                       remesh_global, transfer_state)
from icepack2_tools.geometry import sample_to_geometry  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MESH_DIR = os.path.join(HERE, "..", "mesh")
DATA_DIR = os.path.join(HERE, "..", "data")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("checkpoint")
    ap.add_argument("--out-checkpoint", required=True)
    ap.add_argument("--out-mesh", default=None, help="default: <old basename>_adaptN.msh in antarctica/mesh")
    ap.add_argument("--rebuild-aref", action="store_true")
    ap.add_argument("--no-remesh", action="store_true",
                    help="identity test: keep the old mesh (copied under the new name) and only run the transfer")
    args = ap.parse_args()
    t0 = time.time()
    cfg = AdaptMeshConfig.from_env()
    PETSc.Sys.Print(f"adapt: config {cfg}")

    with fd.CheckpointFile(args.checkpoint, "r") as chk:
        mesh = chk.load_mesh()
        attrs = {k: chk.get_attr("/", k) for k in
                 ("mesh_basename", "buffer_m", "geometry_space", "raster_sample", "adapt_count", "t_yr")
                 if chk.has_attr("/", k)}
        H = chk.load_function(mesh, name="thickness")
        b = chk.load_function(mesh, name="bed")
        try:
            u = chk.load_function(mesh, name="velocity")
        except Exception:
            u = None
    basename = str(attrs.get("mesh_basename", "")).replace(".msh", "")
    if not basename:
        raise RuntimeError("checkpoint has no mesh_basename attribute; cannot find its .msh/sidecar")
    old_msh = os.path.join(MESH_DIR, basename + ".msh")
    old_sidecar = os.path.join(MESH_DIR, f"boundary_ids_{basename}.json")
    for f in (old_msh, old_sidecar):
        if not os.path.exists(f):
            raise FileNotFoundError(f)
    # extract_ice_outline() reads ISMIP7_BUFFER_M: the new mesh must use the
    # old mesh's buffer, not whatever the environment says.
    os.environ["ISMIP7_BUFFER_M"] = str(float(attrs.get("buffer_m", 0.0)))
    k = int(attrs.get("adapt_count", 0)) + 1
    root = basename.split("_adapt")[0]
    out_msh = args.out_mesh or os.path.join(MESH_DIR, f"{root}_adapt{k}.msh")
    new_basename = os.path.splitext(os.path.basename(out_msh))[0]
    PETSc.Sys.Print(f"adapt: {basename} (t={attrs.get('t_yr', '?')}) -> {new_basename}")

    if args.no_remesh:
        import shutil
        if mesh.comm.rank == 0:
            shutil.copy(old_msh, out_msh)
            shutil.copy(old_sidecar, os.path.join(MESH_DIR, f"boundary_ids_{new_basename}.json"))
        mesh.comm.barrier()
        n_ele = mesh.comm.allreduce(FunctionSpace(mesh, "DG", 0).dof_dset.size)
        PETSc.Sys.Print("adapt: --no-remesh, transferring onto a fresh load of the same mesh")
    else:
        h_des, diag = desired_element_size(mesh, cfg, H, b, u=u)
        n_ele = remesh_global(mesh, h_des, cfg, out_msh, old_msh, old_sidecar)
    mesh_new = fd.Mesh(out_msh, name="firedrake_default")

    bm_fn = sorted(glob.glob(os.path.join(DATA_DIR, "bedmachine", "*.nc")))[0]
    method = str(attrs.get("raster_sample", "vertex"))

    def bed_sampler(Q_g, Qc):
        return sample_to_geometry(rasterio.open(f"netcdf:{bm_fn}:bed"), Q_g, Qc, method=method)

    def thickness_sampler(Q_g, Qc):
        return sample_to_geometry(rasterio.open(f"netcdf:{bm_fn}:thickness"), Q_g, Qc,
                                  floor=0.0, method=method)

    audit = transfer_state(args.checkpoint, mesh_new, cfg, args.out_checkpoint,
                           new_basename, bed_sampler, rebuild_aref=args.rebuild_aref,
                           thickness_sampler=thickness_sampler)
    n_old = mesh.comm.allreduce(FunctionSpace(mesh, "DG", 0).dof_dset.size)   # owned cells only
    PETSc.Sys.Print(f"adapt: {n_old} -> {n_ele} cells; bed re-sampled ({method}); "
                    f"volume change {audit['volume_change_pct']:+.4f}%; "
                    f"wrote {args.out_checkpoint} in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
