# -*- coding: utf-8 -*-
"""Render the HUDF template onto a PNG: smooth Sersic vs. clumpy galaxies.

Visual inspection only. Paints the image fields (and clump table) of
``extragalactic.hubble_ultra_deep_field`` onto a common pixel grid using
their ``CRVAL`` offsets, without running an optical train, and saves a
log-stretched two-panel PNG next to this script.

Axes follow the template convention: ``x`` is the offset towards East,
plotted increasing to the right (i.e. mirrored with respect to a
conventional North-up/East-left sky image).

Run from anywhere::

    python docs/notebooks/plot_hubble_ultra_deep_field.py
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import zoom

from scopesim_templates.extragalactic import hubble_ultra_deep_field


HERE = Path(__file__).resolve().parent
PIXEL_SCALE = 0.03       # arcsec / pixel
FOV_ARCSEC = 30.0        # square side
MAG_LIMIT = 28.0
SEED = 42


def _paint(src, npix, pixel_scale):
    canvas = np.zeros((npix, npix))
    centre = (npix - 1) / 2.0
    for field in src.image_fields:
        hdr = field.header
        data = np.asarray(field.data, dtype=float)
        scale = abs(hdr["CDELT1"]) * 3600.0
        if not np.isclose(scale, pixel_scale):
            total = data.sum()
            data = zoom(data, scale / pixel_scale, order=1)
            data *= total / max(data.sum(), 1e-30)
        h, w = data.shape
        x0 = int(round(centre + hdr["CRVAL1"] * 3600.0 / pixel_scale)) - w // 2
        y0 = int(round(centre + hdr["CRVAL2"] * 3600.0 / pixel_scale)) - h // 2
        cx0, cy0 = max(0, x0), max(0, y0)
        cx1, cy1 = min(npix, x0 + w), min(npix, y0 + h)
        if cx0 >= cx1 or cy0 >= cy1:
            continue
        canvas[cy0:cy1, cx0:cx1] += data[cy0 - y0:cy1 - y0, cx0 - x0:cx1 - x0]

    for field in src.table_fields:
        ix = np.round(centre + np.asarray(field.field["x"]) / pixel_scale)
        iy = np.round(centre + np.asarray(field.field["y"]) / pixel_scale)
        ok = (ix >= 0) & (ix < npix) & (iy >= 0) & (iy < npix)
        np.add.at(canvas, (iy[ok].astype(int), ix[ok].astype(int)),
                  np.asarray(field.field["weight"])[ok])
    return canvas


def main(output=HERE / "hubble_ultra_deep_field.png"):
    npix = int(round(FOV_ARCSEC / PIXEL_SCALE))
    kwargs = dict(pixel_scale=PIXEL_SCALE, fov=FOV_ARCSEC, margin=5.0,
                  mag_limit=MAG_LIMIT)
    panels = [
        ("smooth Sersic", hubble_ultra_deep_field(**kwargs)),
        ("with clumps", hubble_ultra_deep_field(clumps=True, seed=SEED,
                                                **kwargs)),
    ]

    images = [_paint(src, npix, PIXEL_SCALE) for _, src in panels]
    vmax = max(im.max() for im in images)
    vmin = vmax * 1e-5
    half = FOV_ARCSEC / 2
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.2), dpi=110)
    for ax, (title, src), im in zip(axes, panels, images):
        ax.imshow(np.clip(im, vmin, None), origin="lower", cmap="afmhot",
                  norm=LogNorm(vmin=vmin, vmax=vmax),
                  extent=[-half, half, -half, half])
        ax.set_title(f"{title}: {src.meta['n_galaxies']} galaxies, "
                     f"{src.meta['n_clumps']} clumps")
        ax.set_xlabel(r"$x$ = $\Delta\alpha\cos\delta$ [arcsec] (East $\rightarrow$)")
    axes[0].set_ylabel(r"$\Delta\delta$ [arcsec] (North $\uparrow$)")
    fig.suptitle(f"HUDF template, m_F160W < {MAG_LIMIT}, "
                 f"{PIXEL_SCALE}\"/pix (model, no PSF)")
    fig.tight_layout()
    fig.savefig(output)
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
