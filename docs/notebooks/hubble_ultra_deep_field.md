---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Hubble Ultra Deep Field

`extragalactic.hubble_ultra_deep_field()` renders the real HUDF galaxy
population (van der Wel et al. 2012 F160W Sersic fits, cross-matched with
Skelton et al. 2014 3D-HST redshifts and masses; 3152 galaxies within
2.5 arcmin of the field centre) as one flux-normalised Sersic image per
galaxy. Optionally, a literature-motivated fraction of the galaxies is
"clumpified": 10-25 per cent of their light is moved into power-law
distributed point-source clumps, conserving the total flux.

## Quick example

```{code-cell} ipython3
from scopesim_templates.extragalactic import hubble_ultra_deep_field

src = hubble_ultra_deep_field(fov=20.0, mag_limit=27.0)
src.meta["n_galaxies"], src.meta["total_flux_F160W_mag"]
```

Each image field is placed with `CRVAL1/2` = offset from
(`ra_center`, `dec_center`) in degrees, with `x` towards East and `y`
towards North; the absolute field centre is stored in `src.meta`. Each
image sums to the galaxy's catalogue flux, $10^{-0.4\,m_{F160W}}$, and
its spectrum is normalised to 0 mag(AB) in `filter_name`.

With clumps (reproducible with `seed`):

```{code-cell} ipython3
src = hubble_ultra_deep_field(fov=16.384, pixel_scale=0.004, clumps=True,
                              seed=42)
src.meta["n_clumpy_galaxies"], src.meta["n_clumps"]
```

## Picking a field

The helpers in `extragalactic.deep_field` expose the catalogue directly:

```{code-cell} ipython3
from astropy import units as u
from astropy.coordinates import SkyCoord
from scopesim_templates.extragalactic.deep_field import (
    HUDF_RA, HUDF_DEC, load_hudf_catalogue, richest_field)

cat = load_hudf_catalogue()
x0, y0, score = richest_field(cat, size=16.384, z_range=(1, 3))
centre = SkyCoord(HUDF_RA, HUDF_DEC, unit="deg").spherical_offsets_by(
    x0 * u.arcsec, y0 * u.arcsec)
src = hubble_ultra_deep_field(fov=16.384, ra_center=centre.ra.deg,
                              dec_center=centre.dec.deg, margin=4.0)
```

## Model image

[`plot_hubble_ultra_deep_field.py`](./plot_hubble_ultra_deep_field.py)
paints the source fields onto a common grid (no optical train):

```{figure} hubble_ultra_deep_field.png
:alt: HUDF template, smooth and clumpy
:width: 700px

A 30 arcsec HUDF patch at 30 mas/pixel, smooth (left) and with clumps
(right). East is to the right (see the module docstring for the
convention).
```

Please cite van der Wel et al. 2012 (ApJS 203, 24) and Skelton et al. 2014
(ApJS 214, 24) when using this template.
