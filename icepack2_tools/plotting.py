"""
Visualization utilities for icepack2 results.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def mesh_plot(mesh, figsize=(20, 14), title=None, save_path=None, zooms=None):
    """Plot mesh with element size colormap and optional zoom panels.

    Parameters
    ----------
    mesh : firedrake.Mesh
        The mesh to plot.
    zooms : list of (name, xlim, ylim) tuples
        Zoom regions to show in additional panels.
    """
    coords = mesh.coordinates.dat.data_ro
    cells = mesh.coordinates.cell_node_map().values
    x, y = coords[:, 0] / 1e3, coords[:, 1] / 1e3
    triang = mtri.Triangulation(x, y, cells)

    # Element sizes
    areas = np.array(
        [
            0.5
            * abs(
                (coords[tri[1], 0] - coords[tri[0], 0])
                * (coords[tri[2], 1] - coords[tri[0], 1])
                - (coords[tri[2], 0] - coords[tri[0], 0])
                * (coords[tri[1], 1] - coords[tri[0], 1])
            )
            for tri in cells
        ]
    )
    h_elem = np.sqrt(areas) / 1e3

    if zooms is None:
        zooms = [
            ("Ross", (-500, 400), (-1400, -500)),
            ("PIG/Thwaites", (-1800, -1350), (-600, -100)),
        ]

    n_panels = 1 + len(zooms) + 1  # continent + zooms + histogram
    ncols = min(3, n_panels)
    nrows = (n_panels + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_2d(axes)

    # Continent view
    ax = axes.flat[0]
    tc = ax.tripcolor(triang, h_elem, cmap="viridis_r", vmin=0, vmax=15, shading="flat")
    ax.set_aspect("equal")
    ax.set_xlim(-2800, 2800)
    ax.set_ylim(-2800, 2800)
    ax.set_title(
        f"Element size\n({mesh.num_vertices()} verts, " f"{mesh.num_cells()} cells)"
    )
    plt.colorbar(tc, ax=ax, shrink=0.6, label="km")

    # Zoom panels
    for i, (name, xlim, ylim) in enumerate(zooms):
        ax = axes.flat[1 + i]
        ax.triplot(triang, color="k", lw=0.1, alpha=0.5)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.set_title(name)

    # Histogram
    ax = axes.flat[-1]
    ax.hist(h_elem, bins=100, color="steelblue", edgecolor="k", lw=0.3)
    ax.set_xlabel("Element size (km)")
    ax.set_ylabel("Count")
    ax.set_title("Size distribution")
    ax.set_xlim(0, 20)

    if title:
        fig.suptitle(title, fontsize=14, y=1.0)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig
