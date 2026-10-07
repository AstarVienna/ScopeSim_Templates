# -*- coding: utf-8 -*-
"""Render the HUDF analogue Source onto a single mosaic PNG.

Visual inspection only — verifies that
``extragalactic.hubble_ultra_deep_field`` builds a sensible scene at
the chosen pixel scale. Each catalogue galaxy contributes one image
HDU with its own WCS; this script paints them onto a common pixel
grid centred on the catalogue centre and saves a log-stretched PNG.

Run from the repo root::

    python scopesim_templates/tests/visual_inspection/plot_hubble_ultra_deep_field.py
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.table import Table
from astropy.wcs import WCS

from scopesim_templates.extragalactic import hubble_ultra_deep_field
from scopesim_templates.extragalactic.deep_field import _DEFAULT_CATALOGUE


PIXEL_SCALE = 0.06       # arcsec / pixel
FOV_ARCSEC = 60.0        # square side
MAG_LIMIT = 27.0
OUT_PATH = Path(__file__).resolve().parent / "hubble_ultra_deep_field.png"


def _paint_field(field, canvas, canvas_wcs, flux):
    """Drop one galaxy field onto the shared canvas via WCS-aware shift.

    ``flux`` is the linear (10^(-0.4*mag)) scaling that compensates for
    the per-stamp sum-normalisation done inside galaxy().
    """
    data = field.data * float(flux)
    h, w = data.shape
    src_wcs = WCS(field.header)
    cx, cy = (w - 1) / 2.0, (h - 1) / 2.0
    sky = src_wcs.pixel_to_world(cx, cy)
    px, py = canvas_wcs.world_to_pixel(sky)
    ix, iy = int(round(float(px))), int(round(float(py)))

    x0 = ix - w // 2
    y0 = iy - h // 2
    x1 = x0 + w
    y1 = y0 + h
    cx0 = max(0, x0)
    cy0 = max(0, y0)
    cx1 = min(canvas.shape[1], x1)
    cy1 = min(canvas.shape[0], y1)
    if cx0 >= cx1 or cy0 >= cy1:
        return

    sx0 = cx0 - x0
    sy0 = cy0 - y0
    sx1 = sx0 + (cx1 - cx0)
    sy1 = sy0 + (cy1 - cy0)
    canvas[cy0:cy1, cx0:cx1] += data[sy0:sy1, sx0:sx1]


def main():
    print(f"Building HUDF Source (fov={FOV_ARCSEC}\", mag_limit={MAG_LIMIT})...")
    src = hubble_ultra_deep_field(
        pixel_scale=PIXEL_SCALE,
        fov=FOV_ARCSEC,
        mag_limit=MAG_LIMIT,
    )
    n_gal = src.meta["n_galaxies"]
    print(f"  {n_gal} galaxies")

    npix = int(round(FOV_ARCSEC / PIXEL_SCALE))
    canvas = np.zeros((npix, npix), dtype=float)

    # Canvas WCS centred on the catalogue centre.
    canvas_wcs = WCS(naxis=2)
    canvas_wcs.wcs.crpix = [(npix + 1) / 2, (npix + 1) / 2]
    canvas_wcs.wcs.cdelt = [-PIXEL_SCALE / 3600.0, PIXEL_SCALE / 3600.0]
    canvas_wcs.wcs.crval = [src.meta["ra_center"], src.meta["dec_center"]]
    canvas_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    canvas_wcs.wcs.cunit = ["deg", "deg"]

    # Pull the same row order back from the catalogue to recover each
    # galaxy's apparent F160W magnitude (Source.fields holds the
    # normalised stamps; the magnitudes live on the catalogue, not on
    # the field headers).
    tbl = Table.read(_DEFAULT_CATALOGUE)
    half = FOV_ARCSEC / 2 / 3600.0
    cra, cdec = src.meta["ra_center"], src.meta["dec_center"]
    cdec_rad = np.deg2rad(cdec)
    sel = (
        (np.abs((np.asarray(tbl["ra"]) - cra) * np.cos(cdec_rad)) <= half)
        & (np.abs(np.asarray(tbl["dec"]) - cdec) <= half)
        & (np.asarray(tbl["m_F160W"]) <= MAG_LIMIT)
    )
    mags = np.asarray(tbl["m_F160W"])[sel]
    if len(mags) != n_gal:
        raise RuntimeError(
            f"catalogue selection ({len(mags)}) != Source.n_galaxies ({n_gal})")
    fluxes = 10.0 ** (-0.4 * mags)

    for i, (field, flux) in enumerate(zip(src.fields, fluxes)):
        _paint_field(field, canvas, canvas_wcs, flux)
        if (i + 1) % 200 == 0:
            print(f"  painted {i + 1}/{n_gal}")

    floor = max(canvas[canvas > 0].min() * 1e-2, 1e-12) if (canvas > 0).any() else 1e-12
    fig, ax = plt.subplots(figsize=(8, 8), dpi=120)
    ax.imshow(
        np.clip(canvas, floor, None),
        origin="lower",
        norm=LogNorm(vmin=floor, vmax=canvas.max() or 1.0),
        cmap="afmhot",
        extent=[-FOV_ARCSEC / 2, FOV_ARCSEC / 2,
                -FOV_ARCSEC / 2, FOV_ARCSEC / 2],
    )
    ax.set_xlabel(r"$\Delta$RA [arcsec]")
    ax.set_ylabel(r"$\Delta$Dec [arcsec]")
    ax.set_title(
        f"HUDF analogue: {n_gal} galaxies, m_F160W < {MAG_LIMIT}\n"
        f"pixel scale = {PIXEL_SCALE}\"/pix, FOV = {FOV_ARCSEC}\""
    )
    fig.tight_layout()
    fig.savefig(OUT_PATH)
    print(f"Saved {OUT_PATH}")


if __name__ == "__main__":
    main()
