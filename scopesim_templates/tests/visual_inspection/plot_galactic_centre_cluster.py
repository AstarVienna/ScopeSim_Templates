# -*- coding: utf-8 -*-
"""Render a 20-year GIF animation of stellar.galactic_centre().

Each frame is built from a fresh scopesim.Source produced by
``galactic_centre(time)`` and reads positions, RVs and weights from
the Source.fields table.

Run from anywhere::

    python scopesim_templates/tests/visual_inspection/plot_galactic_centre_cluster.py

The output GIF is written next to this script.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from astropy import units as u
from astropy.time import Time

from scopesim_templates.stellar import galactic_centre


HERE = Path(__file__).resolve().parent


def main(output: str = str(HERE / "galactic_centre_evolution.gif"),
         start: str = "2024-01-01",
         duration_years: float = 20.0,
         n_frames: int = 41,
         fps: int = 8,
         rv_clip_kms: float = 1500.0,
         axis_limit_arcsec: float = 1.0):
    t0 = Time(start)
    times = t0 + np.linspace(0.0, duration_years, n_frames) * u.yr

    def source_at(t):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return galactic_centre(t)

    src0 = source_at(times[0])
    f0 = src0.fields[0].field
    weight = np.asarray(f0["weight"])
    # Sizes scaled by brightness; clip extremes for legibility
    sizes = np.clip(120.0 * weight / weight.max(), 4.0, 200.0)

    fig, ax = plt.subplots(figsize=(7.0, 7.0))
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    ax.set_xlim(axis_limit_arcsec, -axis_limit_arcsec)  # east-left
    ax.set_ylim(-axis_limit_arcsec, axis_limit_arcsec)
    ax.set_xlabel(r"$\Delta\alpha\cos\delta$ [arcsec]  (East $\leftarrow$)")
    ax.set_ylabel(r"$\Delta\delta$ [arcsec]  ($\uparrow$ North)")
    ax.grid(alpha=0.2, linestyle=":")
    ax.scatter([0], [0], marker="+", s=200, c="yellow",
               linewidths=2, zorder=5, label="Sgr A*")

    sc = ax.scatter(
        np.asarray(f0["x"]), np.asarray(f0["y"]),
        s=sizes, c=np.asarray(f0["rv"]),
        cmap="RdBu_r", vmin=-rv_clip_kms, vmax=rv_clip_kms,
        edgecolors="white", linewidths=0.3, zorder=4,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Radial velocity [km/s]")

    date_text = ax.text(
        0.02, 0.98, "", transform=ax.transAxes,
        va="top", ha="left", color="white",
        fontsize=12, family="monospace",
    )
    ax.set_title("Galactic Centre — Gillessen+2017 orbits")
    ax.legend(loc="lower right", facecolor="black",
              edgecolor="white", labelcolor="white")

    def update(idx):
        src = source_at(times[idx])
        f = src.fields[0].field
        sc.set_offsets(np.column_stack([np.asarray(f["x"]),
                                        np.asarray(f["y"])]))
        sc.set_array(np.asarray(f["rv"]))
        date_text.set_text(times[idx].iso[:10])
        return sc, date_text

    anim = FuncAnimation(fig, update, frames=n_frames,
                         interval=1000.0 / fps, blit=False)
    anim.save(output, writer=PillowWriter(fps=fps), dpi=110)
    plt.close(fig)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
