# -*- coding: utf-8 -*-
"""Render a 20-year GIF animation of stellar.globular_cluster().

Each frame is built from a fresh scopesim.Source produced by
``globular_cluster(..., time=t, seed=fixed)`` and reads positions, RVs
and weights from the Source.fields table.

Run from anywhere::

    python docs/notebooks/plot_globular_cluster.py

The output GIF is written next to this script, where the companion
notebook ``globular_cluster.md`` picks it up.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

from scopesim_templates.stellar import globular_cluster


HERE = Path(__file__).resolve().parent


def main(output: str = str(HERE / "globular_cluster_evolution.gif"),
         density: float = 750.0,
         fov: float = 2.0,
         distance_modulus: float = 13.6,
         imbh_mass: float = 1e5,
         seed: int = 42,
         duration_years: float = 20.0,
         n_frames: int = 41,
         fps: int = 8,
         rv_clip_kms: float = 500.0,
         axis_limit_arcsec: float = 1.0):

    times_yr = np.linspace(0.0, duration_years, n_frames)

    def source_at(t):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return globular_cluster(
                density=density, fov=fov,
                distance_modulus=distance_modulus,
                imbh_mass=imbh_mass, seed=seed, time=float(t),
            )

    src0 = source_at(0.0)
    f0 = src0.fields[0].field
    weight = np.asarray(f0["weight"])
    # Brightness-scaled marker sizes; clip extremes for legibility
    sizes = np.clip(15.0 * weight / weight.max(), 1.5, 30.0)

    half = float(axis_limit_arcsec)
    fig, ax = plt.subplots(figsize=(7.5, 7.0))
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    ax.set_xlim(half, -half)
    ax.set_ylim(-half, half)
    ax.set_xlabel(r"$\Delta\alpha\cos\delta$ [arcsec]  (East $\leftarrow$)")
    ax.set_ylabel(r"$\Delta\delta$ [arcsec]  ($\uparrow$ North)")
    ax.grid(alpha=0.2, linestyle=":")

    # FOV outline
    fov_half = fov / 2.0
    ax.plot([fov_half, fov_half, -fov_half, -fov_half, fov_half],
            [-fov_half, fov_half, fov_half, -fov_half, -fov_half],
            color="white", linewidth=0.6, alpha=0.4, linestyle="--")

    ax.scatter([0], [0], marker="*", s=200, c="yellow",
               linewidths=1, zorder=5, label="IMBH")

    sc = ax.scatter(
        np.asarray(f0["x"]), np.asarray(f0["y"]),
        s=sizes, c=np.asarray(f0["rv"]),
        cmap="RdBu_r", vmin=-rv_clip_kms, vmax=rv_clip_kms,
        edgecolors="none", zorder=4,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Radial velocity [km/s]")

    elapsed_text = ax.text(
        0.02, 0.98, "", transform=ax.transAxes,
        va="top", ha="left", color="white",
        fontsize=12, family="monospace",
    )
    ax.set_title(
        f"Globular cluster + IMBH ({imbh_mass:.0e} M$_\\odot$, "
        f"DM={distance_modulus:.1f}, seed={seed})"
    )
    ax.legend(loc="lower right", facecolor="black",
              edgecolor="white", labelcolor="white")

    def update(idx):
        src = source_at(times_yr[idx])
        f = src.fields[0].field
        sc.set_offsets(np.column_stack([np.asarray(f["x"]),
                                        np.asarray(f["y"])]))
        sc.set_array(np.asarray(f["rv"]))
        elapsed_text.set_text(f"t = {times_yr[idx]:5.2f} yr")
        return sc, elapsed_text

    anim = FuncAnimation(fig, update, frames=n_frames,
                         interval=1000.0 / fps, blit=False)
    anim.save(output, writer=PillowWriter(fps=fps), dpi=110)
    plt.close(fig)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
