#!/usr/bin/python3
# Nicholas Rathmann and Daniel Shapero, 2024

r"""
Test firedrake interface for specfab.

Assumes a time-constant, non-uniform shear flow in vertical cross-section (xz) domain.
"""
import os
import numpy as np

import rasterio

import firedrake as fd
import icepack

"""
Fabric problem setup
"""

def main():
    if os.path.exists("/N/slate/dlilien"):
        proj_dir = "/N/slate/dlilien"
    elif os.path.exists("/Volumes/slate"):
        proj_dir = "/Volumes/slate"
    else:
        raise FileNotFoundError("Mount slate!")

    smb_fn = os.path.join(proj_dir, "greenland_general/climate/GRN11_RACMO24p1_data/smbgl_monthlyS_GRN11_RACMO24p1_200601_201512_mean.tif")
    bed_fn = "netcdf:" + os.path.join(proj_dir, "greenland_general/bedmachine/BedMachineGreenland-v5.nc:bed")
    thick_fn = "netcdf:" + os.path.join(proj_dir, "greenland_general/bedmachine/BedMachineGreenland-v5.nc:thickness")
    vx_fn = os.path.join(proj_dir, "greenland_general/velocity/multiyear/greenland_vel_mosaic250_vx_v1.tif")
    vy_fn = os.path.join(proj_dir, "greenland_general/velocity/multiyear/greenland_vel_mosaic250_vy_v1.tif")
    ex_fn = os.path.join(proj_dir, "greenland_general/velocity/multiyear/greenland_vel_mosaic250_ex_v1.tif")
    ey_fn = os.path.join(proj_dir, "greenland_general/velocity/multiyear/greenland_vel_mosaic250_ey_v1.tif")
    vx = rasterio.open(vx_fn, masked=True)
    vy = rasterio.open(vy_fn, masked=True)

    ex = rasterio.open(ex_fn, masked=True)
    ey = rasterio.open(ey_fn, masked=True)

    smb_arr = rasterio.open(smb_fn)
    thick_arr = rasterio.open(thick_fn)
    bed_arr = rasterio.open(bed_fn)

    for mesh_dir in ["promice_meshes", "simple_meshes", "meshes"]:
        mesh_fn_template = mesh_dir + "/greenland_{:d}_{:d}.msh"
        out_fn_template = mesh_dir + "/greenland_{:d}_{:d}.h5"

        for lc in [250 * 2 ** i for i in range(8)]:
            print("Doing {:d}...".format(lc))
            mesh_fn = mesh_fn_template.format(10 * lc, lc)
            mesh2d = fd.Mesh(mesh_fn)
            mesh = fd.ExtrudedMesh(mesh2d, layers=1)
            mesh.name = "greenland"

            V = fd.VectorFunctionSpace(mesh, "CG", 1, dim=2, vfamily="R", vdegree=0)
            Q = fd.FunctionSpace(mesh, "CG", 1, vfamily="R", vdegree=0) # for projecting scalar fabric measures

            smb = fd.Function(Q).interpolate(fd.max_value(icepack.interpolate(smb_arr, Q), fd.Constant(0.0)))  # in kg m-2
            u = icepack.interpolate((vx, vy), V, rplc_ndv=[(-2.0e9, 0.0), (-2.0e8, 0.0), (2.0e8, 0.0)], fillvalue=0.0)
            sigma_u = icepack.interpolate((ex, ey), V, rplc_ndv=[(-2.0e9, 1.0e4), (-2.0e8, 1.0e4), (2.0e8, 1.0e4)], fillvalue=1.04)
            H = fd.Function(Q).interpolate(fd.max_value(icepack.interpolate(thick_arr, Q), fd.Constant(10.0)))
            b = icepack.interpolate(bed_arr, Q)

            print("Check on the minimum velocity to ensure that we masked:", np.min(u.dat.data[:]))

            out_fn = out_fn_template.format(10 * lc, lc)
            with fd.CheckpointFile(out_fn, "w") as chk:
                chk.save_mesh(mesh)
                chk.save_function(u, name="u")
                chk.save_function(sigma_u, name="sigma_u")
                chk.save_function(smb, name="smb")
                chk.save_function(H, name="H")
                chk.save_function(b, name="b")

if __name__ == "__main__":
    main()
