#!/usr/bin/env python
# coding: utf-8
import os
import firedrake
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import sys
import xarray as xr
import rioxarray as rxr
from shapely.geometry import LineString
from shapely import simplify

# velocity_filename = '../../QGreenland_v3.0.0/Glaciology/Ice sheet velocity/GrIMP/Annual ice sheet velocity magnitude 2021 (200m)/grimp_annual_vv_2021.tif'
vx_fn = "/Volumes/LaCie/Data/greenland_general/velocity/multiyear/greenland_vel_mosaic250_vx_v1.tif"
vy_fn = "/Volumes/LaCie/Data/greenland_general/velocity/multiyear/greenland_vel_mosaic250_vy_v1.tif"
vx = rxr.open_rasterio(vx_fn, masked=True)
vy = rxr.open_rasterio(vy_fn, masked=True)
#vx = vx_f.VX
#vy = vy_f.VY

dudx = vx.differentiate("x")
dudy = vx.differentiate("y")
dvdx = vy.differentiate("x")
dvdy = vy.differentiate("y")

mag_eps = np.abs(dudx) + np.abs(dvdy) + 0.5 * (np.abs(dudy) + np.abs(dvdx))
plt.imshow(mag_eps[0, :, :], vmin=0, vmax=1, cmap="PuOr", extent=[vx.x[0], vx.x[-1], vx.y[-1], vx.y[0]])
plt.colorbar()

plt.figure()


mag_eps = (mag_eps * 20 + np.sqrt(vx ** 2.0 + vy ** 2.0) / 2000) ** 0.75
plt.imshow(mag_eps[0, :, :], vmin=0, vmax=1, cmap="PuOr", extent=[vy.x[0], vy.x[-1], vy.y[-1], vy.y[0]])
plt.colorbar()


outline0 = gpd.read_file("/Volumes/LaCie/Data/greenland_general/drainage_basins/promice_outline/02-PROMICE-2022-IceMask-polygon-v3.gpkg").to_crs("EPSG:3413")
# fig, ax = plt.subplots()
# outline0.plot(ax=ax)
outline0.geometry[0] = simplify(outline0.geometry[0], 250)
outline0.to_file("test_250.gpkg", driver="GPKG")
# outline0.plot(ax=ax, color="C1")
# plt.show()

outline1 = gpd.read_file("/Volumes/LaCie/Data/greenland_general/drainage_basins/greenland_ice_sheet_comb_single_simp.gpkg").to_crs("EPSG:3413")
outline2 = gpd.read_file("simple_meshes/simple_polygon_of_greenland.gpkg").to_crs("EPSG:3413")
outlines = [outline0, outline1, outline2]


meps = mag_eps[0, ::2, ::2]
meps = meps.fillna(1.0e-8)
fig, ax = plt.subplots()
ax.imshow(meps, vmin=0, vmax=1, cmap="PuOr", extent=[meps.x[0], meps.x[-1], meps.y[-1], meps.y[0]])

for outline, out_dir, stride in zip(outlines, ["promice_meshes", "meshes", "simple_meshes"], [1, 2, 32]):
    coord_arr_rough = np.array(outline.geometry[0].exterior.xy).T[:, :]

    # Dork around with refining edges because we could be really rough
    coord_arr = np.zeros(((coord_arr_rough.shape[0] - 1) * stride + 1, 2))
    coord_arr[0::stride, :] = coord_arr_rough
    for i in range(1, stride):
        coord_arr[i::stride, :] = (coord_arr_rough[:-1, :] * (stride - i) + coord_arr_rough[1:, :] * i) / stride

    x = coord_arr[:, 0]
    y = coord_arr[:, 1]
    eval_x = xr.DataArray((x[1:] + x[:-1]) / 2.0 , dims="v")
    eval_y = xr.DataArray((y[1:] + y[:-1]) / 2.0 , dims="v")

    normal_x = -(y[1:] - y[:-1])
    normal_y = (x[1:] - x[:-1])
    norm = np.sqrt(normal_x ** 2.0 + normal_y ** 2.0)
    normal_x /= norm
    normal_y /= norm

    vx_out = vx.interp(x=eval_x, y=eval_y)
    vy_out = vy.interp(x=eval_x, y=eval_y)

    normal_vel = vx_out * normal_x + vy_out * normal_y
    nv = normal_vel.values.flatten()
    nv[np.isnan(nv)] = 9999

    fig, ax = plt.subplots()
    ax.imshow(np.sqrt(vx.values**2.0 + vy.values ** 2.0)[0, :, :], vmin=0, vmax=25, extent=[meps.x[0], meps.x[-1], meps.y[-1], meps.y[0]])
    ax.plot(x, y, color='k')
    ax.scatter(eval_x.values, eval_y.values, c=nv, vmin=0, vmax=25)

    fig, ax = plt.subplots()
    ax.imshow(np.sqrt(vx.values**2.0 + vy.values ** 2.0)[0, :, :], vmin=0, vmax=25, extent=[meps.x[0], meps.x[-1], meps.y[-1], meps.y[0]])
    ax.plot(x, y, color='gray')
    ax.quiver(eval_x.values, eval_y.values, normal_x, normal_y, scale=50.0)
    # ax.quiver(eval_x.values, eval_y.values, vx_out, vy_out, scale=1.0, color="w")

    cutoff_vel = 10.0
    this_line_is_outlet = nv[0] > cutoff_vel
    coords = []
    names = []
    outline_lc_map = []
    start_ind = 0
    for i in range(coord_arr.shape[0] - 1):
        if ((nv[i] > cutoff_vel) != this_line_is_outlet) or (i == (coord_arr.shape[0] - 2)):
            if this_line_is_outlet:
                names.append("Outlet {:d}".format(len(coords)))
            else:
                names.append("Other {:d}".format(len(coords)))
            coords.append(coord_arr[start_ind:i, :])
            outline_lc_map.append(this_line_is_outlet)
            start_ind = i
            this_line_is_outlet = nv[i] > cutoff_vel

    df = gpd.GeoDataFrame(data={"Name": names}, geometry=[LineString(np.vstack((coords[i - 1][-1, :], coords[i]))) for i in range(len(coords))], crs="EPSG:3413")
    df_fn = "boundaries_greenland.gpkg"
    if not os.path.exists(df_fn):
        df.to_file(df_fn, driver="GPKG")

    fig, ax = plt.subplots()
    ax.imshow(np.sqrt(vx.values**2.0 + vy.values ** 2.0)[0, :, :], vmin=0, vmax=25, extent=[meps.x[0], meps.x[-1], meps.y[-1], meps.y[0]])
    for name, c in zip(names, coords):
        color = "k"
        if name == "Other":
            color = "gray"
        ax.plot(c[:, 0], c[:, 1], color=color)
    ax.scatter(eval_x.values, eval_y.values, c=nv, vmin=0, vmax=25, zorder=999999)

    # names = ["Border"]
    # coords = [np.array(outline.geometry[0].exterior.xy).T[:-1, :]]

    cuts = [np.array(geom.xy).T[:-1, :] for geom in outline.geometry[0].interiors]
    cuts = [cut for cut in cuts if cut.shape[0] > 3]
    # cuts = cuts[:2]
    cuts_lc = ["Fine" for i in cuts]


    # cuts = []
    # cuts_lc = []



    def make_mesh(rough_targ, fine_targ, outline_lc):
        fig, ax = plt.subplots()
        ax.plot(x, y, color='gray')
        outline_list = []
        cut_list = []


        all_pts = []

        pt_num = 1
        line_num = 1
        for line_index, line_coords in enumerate(coords):
            outline_list.append([])
            current_outline = outline_list[-1]
            first_pt = pt_num
            for pt in line_coords:
                gmsh.model.geo.addPoint(pt[0], pt[1], 0, outline_lc[line_index], pt_num)
                all_pts.append([pt[0], pt[1]])
                pt_num += 1
            if line_index == 0:
                if line_index == len(coords) - 1:
                    print("Just a loop")
                    line = np.hstack((np.arange(first_pt, pt_num), [1]))
                else:
                    line = np.arange(first_pt, pt_num)
            elif line_index == len(coords) - 1:
                line = np.hstack(([first_pt - 1], *np.arange(first_pt, pt_num), [1]))
            else:
                line = np.hstack(([first_pt - 1], *np.arange(first_pt, pt_num)))
            for i in range(len(line) - 1):
                gmsh.model.geo.addLine(line[i], line[i + 1], line_num)
                current_outline.append(line_num)
                line_num += 1

            ax.plot(np.array(all_pts)[np.array(line) - 1, 0], np.array(all_pts)[np.array(line) - 1, 1])
        last_outline = line_num - 1

        for line_index, line_coords in enumerate(cuts):
            cut_list.append([])
            current_cut = cut_list[-1]
            first_pt = pt_num
            if cuts_lc[line_index] == "Fine":
                lc = fine_targ
            else:
                lc = rough_targ

            for pt in line_coords:
                gmsh.model.geo.addPoint(pt[0], pt[1], 0, lc, pt_num)
                pt_num += 1
            line = np.hstack((np.arange(first_pt, pt_num), first_pt))
            for i in range(len(line) - 1):
                gmsh.model.geo.addLine(line[i], line[i + 1], line_num)
                current_cut.append(line_num)
                line_num += 1

        last_cut = line_num - 1

        outline_num = last_cut + 1
        cut_nums = np.arange(outline_num + 1, outline_num + len(cuts) + 1)
        plane_num = last_cut + 1

        gmsh.model.geo.addCurveLoop(np.arange(1, last_outline + 1), outline_num)
        for i, current_cut in enumerate(cut_list):
            gmsh.model.geo.addCurveLoop(current_cut, cut_nums[i])

        gmsh.model.geo.addPlaneSurface([outline_num] + [i for i in cut_nums], plane_num)

        for i, lines in enumerate(outline_list):
            gmsh.model.geo.addPhysicalGroup(1, lines, name=names[i - 1])
        for i, lines in enumerate(cut_list):
            gmsh.model.geo.addPhysicalGroup(1, lines, name="Cut {:d}".format(i))

        gmsh.model.geo.addPhysicalGroup(2, [plane_num], name="Surf")

        gmsh.model.geo.synchronize()
        def meshSizeCallback(dim, tag, x, y, z, lc):
            mep = meps.interp(x=x, y=y, method='nearest').values
            return max(min(rough_targ, fine_targ / mep), fine_targ)

        # gmsh.model.mesh.setSizeCallback(meshSizeCallback)

        # gmsh.model.mesh.generate(1)
        gmsh.model.mesh.generate(2)


    class Mesh:
        def __init__(self):
            self.vtags, vxyz, _ = gmsh.model.mesh.getNodes()
            self.vxyz = vxyz.reshape((-1, 3))
            vmap = dict({j: i for i, j in enumerate(self.vtags)})
            self.triangles_tags, evtags = gmsh.model.mesh.getElementsByType(2)
            evid = np.array([vmap[j] for j in evtags])
            self.triangles = evid.reshape((self.triangles_tags.shape[-1], -1))

    def compute_size_field(nodes, triangles, rough_targ, fine_targ):
        vxyz = nodes[triangles].mean(axis=1)
        x = xr.DataArray(vxyz[:, 0], dims="H")
        y = xr.DataArray(vxyz[:, 1], dims="H")

        tsize = fine_targ / meps.interp(x=x, y=y).values.flatten()
        tsize[tsize < fine_targ] = fine_targ
        tsize[tsize > rough_targ] = rough_targ
        return tsize


    fn_template = out_dir + "/greenland_{:d}_{:d}"

    for lc in [250 * 2 ** i for i in range(8)]:
        rough_targ, fine_targ = lc * 10, lc

        gmsh.initialize(sys.argv)

        outline_lc = [fine_targ if olc else rough_targ for olc in outline_lc_map]

        fn_base = fn_template.format(rough_targ, fine_targ)

        gmsh.model.add(fn_base + "_orig")
        make_mesh(rough_targ, fine_targ, outline_lc)
        gmsh.write(fn_base + "_raw.msh")
        mesh = Mesh()

        sf_ele = compute_size_field(mesh.vxyz, mesh.triangles, rough_targ, fine_targ)
        sf_view = gmsh.view.add("mesh size field")
        gmsh.view.addModelData(sf_view, 0, fn_base + "_orig", "ElementData", mesh.triangles_tags, sf_ele[:, None])
        gmsh.view.write(sf_view, fn_base + "_sf.pos")

        gmsh.model.add(fn_base)
        make_mesh(rough_targ, fine_targ, outline_lc)

        bg_field = gmsh.model.mesh.field.add("PostView")
        gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
        gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

        # In order to compute the mesh sizes from the background mesh only, and
        # disregard any other size constraints, one can set:
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
        gmsh.write(fn_base + ".msh")
        gmsh.finalize()

        os.remove(fn_base + "_raw.msh")
        os.remove(fn_base + "_sf.pos")


    if True:
        fn_base = fn_template.format(20000, 2000)
        mesh2d = firedrake.Mesh(fn_base + ".msh")
        fig, ax = plt.subplots(figsize=(14, 14))
        firedrake.triplot(mesh2d, axes=ax)
        ax.axis("equal")
        ax.legend()
        fig.savefig(fn_base + ".png", dpi=300)
