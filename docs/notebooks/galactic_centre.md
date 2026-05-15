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

# Galactic Centre

`stellar.galactic_centre()` builds a `scopesim.Source` of the Sgr A*
S-star field at any user-supplied epoch. Stellar orbital elements come
from the vendored Gillessen et al. 2017 catalogue
([CDS J/ApJ/837/30](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/837/30));
positions are projected onto the sky at `time` via Kepler's equations
(see API docs). Each star gets one of two spectral types — `"B0V"` for
early-type stars and `"M8III"` for late-type stars by default — scaled
to its catalogue $K$-band magnitude in `Paranal/NACO.Ks` and shifted to
its line-of-sight velocity at the requested epoch.

## Quick example

```{code-cell} ipython3
from astropy.time import Time
from scopesim_templates.stellar import galactic_centre

src = galactic_centre(Time("2024-01-01"))
src.fields[0].field[:5]
```

The columns `x`, `y` are arcsec offsets from Sgr A*, `rv` is the
line-of-sight velocity in km/s at the requested epoch, `weight` is
$10^{-0.4 K_\mathrm{mag}}$, and `ref` indexes into `src.spectra`.

## 20-year evolution

Running the companion script
[`plot_galactic_centre.py`](./plot_galactic_centre.py) builds one
`Source` per frame (`galactic_centre(time)`) and animates the
projected positions colour-coded by RV. S2 swings visibly past Sgr A*
during the window:

```{figure} galactic_centre_evolution.gif
:alt: 20-year animation of the Gillessen S-star field around Sgr A*
:width: 600px

20-year evolution of the Galactic Centre S-star field. Marker size
encodes $K$-band magnitude; colour encodes line-of-sight velocity at
each frame.
```

To regenerate the GIF yourself:

```bash
python docs/notebooks/plot_galactic_centre.py
```
