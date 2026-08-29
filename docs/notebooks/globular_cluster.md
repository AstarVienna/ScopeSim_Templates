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

# Globular Cluster + IMBH

`stellar.globular_cluster()` builds a fictional Omega-Cen-like cluster
with a central intermediate-mass black hole. Stars are drawn from a
Salpeter IMF; the mamajek mass-to-spectral-type / mass-to-Mv tables
provide the per-star spectrum and absolute V magnitude; each star is
placed on a randomly drawn closed Keplerian orbit around the IMBH so
that its sky position and line-of-sight velocity at the user-supplied
`time` fall out together.

## Quick example

```{code-cell} ipython3
from scopesim_templates.stellar import globular_cluster

src = globular_cluster(density=0.5, fov=20.0, distance_modulus=13.6,
                       imbh_mass=1e4, seed=42)
src.fields[0].field[:5]
```

Columns: `x`, `y` are projected sky offsets in arcsec; `rv` is the
line-of-sight velocity in km/s at `time`; `weight` is
$10^{-0.4\,(M_V + \mu)}$; `ref` indexes into `src.spectra`.

Calling `globular_cluster()` again with the same `seed` and a different
`time` advances every star along its orbit, so successive Sources form
a coherent time series suitable for animation.

## 20-year evolution

The companion script
[`plot_globular_cluster.py`](./plot_globular_cluster.py) builds one
`Source` per frame and animates a 3000-star cluster around a
$10^5\,M_\odot$ IMBH:

```{figure} globular_cluster_evolution.gif
:alt: 20-year animation of a globular cluster with a central IMBH
:width: 600px

20-year evolution of a 3000-star cluster around a $10^5\,M_\odot$
IMBH. Inner-orbit stars complete visible loops over the 20-year
window; outer-orbit stars barely move. Colour encodes line-of-sight
velocity at each frame.
```

To regenerate the GIF yourself:

```bash
python docs/notebooks/plot_globular_cluster.py
```
