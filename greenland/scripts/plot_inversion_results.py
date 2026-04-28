#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import firedrake
import matplotlib.pyplot as plt
import numpy as np
from icepackaccs import extract_surface


mesh_dir = "meshes"
fn_template = mesh_dir + "/greenland_{:d}_{:d}.h5"

if mesh_dir == "meshes":
    name = "detailed"
else:
    name = "simple"

all_lcs = np.array([250 * 2 ** i for i in range(0, 8)])

out_template = "inv_{:s}_{:d}.h5"

for lc in all_lcs:
    cache_fn = os.path.join("cached_results", out_template.format(name, lc))
    if not os.path.exists(cache_fn):
        continue
    fn = fn_template.format(lc * 10, lc)
    with firedrake.CheckpointFile(fn, "r") as chk:
        mesh = chk.load_mesh(name="greenland")
        u_obs = chk.load_function(mesh, name="u")
        sigma_u = chk.load_function(mesh, name="sigma_u")
        smb = chk.load_function(mesh, name="smb")
        H_in = chk.load_function(mesh, name="H")
        b_in = chk.load_function(mesh, name="b")


    with firedrake.CheckpointFile(cache_fn, "r") as chk:
        mesh = chk.load_mesh("greenland")
        C = chk.load_function(mesh, name="C")
        u_opt = chk.load_function(mesh, name="u_opt")
        u_obs = chk.load_function(mesh, name="u_obs")

    Q = firedrake.FunctionSpace(
        mesh, "CG", 2, vfamily="R", vdegree=0
    )

    fig, ax = plt.subplots(1, 3, figsize=(14, 7))
    colors = firedrake.tripcolor(extract_surface(firedrake.Function(Q).interpolate(firedrake.sqrt(u_obs[0] ** 2.0 + u_obs[1] ** 2.0))), vmin=0, vmax=1000, axes=ax[0])
    plt.colorbar(colors, ax=ax[0], extend='max', label="Flow speed [m/yr]")
    colors = firedrake.tripcolor(extract_surface(firedrake.Function(Q).interpolate(firedrake.sqrt(u_opt[0] ** 2.0 + u_opt[1] ** 2.0))), vmin=0, vmax=1000, axes=ax[1])
    plt.colorbar(colors, ax=ax[1], extend='max', label="Flow speed [m/yr]")
    ax[0].axis("equal")
    ax[1].axis("equal")
    ax[2].axis("equal")
    colors = firedrake.tripcolor(extract_surface(firedrake.Function(Q).interpolate(firedrake.sqrt(u_obs[0] ** 2.0 + u_obs[1] ** 2.0) - firedrake.sqrt(u_opt[0] ** 2.0 + u_opt[1] ** 2.0))), vmin=-500, vmax=500, axes=ax[2], cmap="PuOr")
    plt.colorbar(colors, ax=ax[2], extend='max', label="Observed - Modeled speed [m/yr]")
    ax[0].set_title("Observed")
    ax[1].set_title("Modeled")
    ax[2].set_title("Difference")
    fig.savefig(os.path.join("figs", "inv_vel_{:s}_{:d}.png".format(name, lc)), dpi=300)
