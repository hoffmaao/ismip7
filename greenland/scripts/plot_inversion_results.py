#!/usr/bin/env python3
"""Plot inversion velocity results for initialized Greenland meshes."""

import argparse
import sys
from pathlib import Path

import numpy as np


DEFAULT_MESH_DIR = Path(__file__).resolve().parent.parent / "meshes"
DEFAULT_LOG_DIR = Path("test_logs")
DEFAULT_FIG_DIR = Path(__file__).resolve().parent.parent / "figs"
MESH_KINDS = ("detailed", "promice", "simple", "buffered")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot inversion velocity results.")
    parser.add_argument(
        "--mesh-dir",
        type=Path,
        default=DEFAULT_MESH_DIR,
        help="Directory containing initialized h5 files.",
    )
    parser.add_argument(
        "--mesh-kind",
        choices=MESH_KINDS,
        default="detailed",
        help="Mesh type to use: detailed, promice, simple, or buffered.",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help="Directory containing cached inversion results.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=DEFAULT_FIG_DIR,
        help="Directory for inversion result figures.",
    )
    return parser.parse_args()


def checkpoint_path(mesh_dir: Path, mesh_kind: str, lc: int) -> Path:
    """Build the checkpoint path matching mesh_greenland.py names."""
    return mesh_dir / f"greenland_{mesh_kind}_{10 * lc}_{lc}.h5"


def plot_results(mesh_dir: Path, mesh_kind: str, log_dir: Path, fig_dir: Path) -> None:
    """Plot observed, modeled, and residual speeds for one mesh kind."""
    import firedrake
    import matplotlib.pyplot as plt
    from icepackaccs import extract_surface

    all_lcs = np.array([250 * 2**i for i in range(0, 8)])
    cache_dir = log_dir / "cached_results"
    out_template = f"inv_{mesh_kind}_{{:d}}.h5"
    fig_dir.mkdir(parents=True, exist_ok=True)

    for lc in all_lcs:
        lc = int(lc)
        cache_fn = cache_dir / out_template.format(lc)
        if not cache_fn.exists():
            continue

        with firedrake.CheckpointFile(str(checkpoint_path(mesh_dir, mesh_kind, lc)), "r") as chk:
            mesh = chk.load_mesh(name="greenland")
            chk.load_function(mesh, name="u")
            chk.load_function(mesh, name="sigma_u")
            chk.load_function(mesh, name="smb")
            chk.load_function(mesh, name="H")
            chk.load_function(mesh, name="b")

        with firedrake.CheckpointFile(str(cache_fn), "r") as chk:
            mesh = chk.load_mesh("greenland")
            u_opt = chk.load_function(mesh, name="u_opt")
            u_obs = chk.load_function(mesh, name="u_obs")

        Q = firedrake.FunctionSpace(mesh, "CG", 2, vfamily="R", vdegree=0)
        obs_speed = firedrake.Function(Q).interpolate(
            firedrake.sqrt(u_obs[0] ** 2.0 + u_obs[1] ** 2.0)
        )
        opt_speed = firedrake.Function(Q).interpolate(
            firedrake.sqrt(u_opt[0] ** 2.0 + u_opt[1] ** 2.0)
        )
        speed_diff = firedrake.Function(Q).interpolate(obs_speed - opt_speed)

        fig, ax = plt.subplots(1, 3, figsize=(14, 7))
        colors = firedrake.tripcolor(
            extract_surface(obs_speed), vmin=0, vmax=1000, axes=ax[0]
        )
        plt.colorbar(colors, ax=ax[0], extend="max", label="Flow speed [m/yr]")
        colors = firedrake.tripcolor(
            extract_surface(opt_speed), vmin=0, vmax=1000, axes=ax[1]
        )
        plt.colorbar(colors, ax=ax[1], extend="max", label="Flow speed [m/yr]")
        colors = firedrake.tripcolor(
            extract_surface(speed_diff), vmin=-500, vmax=500, axes=ax[2], cmap="PuOr"
        )
        plt.colorbar(
            colors, ax=ax[2], extend="max", label="Observed - Modeled speed [m/yr]"
        )

        for axes in ax:
            axes.axis("equal")
        ax[0].set_title("Observed")
        ax[1].set_title("Modeled")
        ax[2].set_title("Difference")
        fig.savefig(fig_dir / f"inv_vel_{mesh_kind}_{lc}.png", dpi=300)
        plt.close(fig)


def main() -> None:
    """Plot inversion results from command-line arguments."""
    args = parse_args()
    sys.argv = [sys.argv[0]]
    plot_results(
        args.mesh_dir.expanduser(),
        args.mesh_kind,
        args.log_dir.expanduser(),
        args.fig_dir.expanduser(),
    )


if __name__ == "__main__":
    main()
