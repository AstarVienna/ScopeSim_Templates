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

# Stellar Module

This module include general functions to work with stars

## Stellar cluster

In the following example we generate a young star cluster with a core radius `r_c=1 pc`, `M=1000 Msun` and located in the LMC (`d=50kpc`)

```{code-cell} ipython3
import numpy as np
from scopesim_templates.stellar import cluster

src = cluster(mass=1E3, distance=50000, core_radius=1)
```

Lets have a look inside the object:

```{code-cell} ipython3
src.fields[0]
```

Here we can see the spatial information is in the form of an `astropy.Table`.

The columns `x` and `y` show the position of each star in `arcsec` relative to the centre of the field of view.

The column `ref` connects each star in this table to a spectrum in the following list:

```{code-cell} ipython3
src.spectra
```

When ScopeSim ingests this Source object, it will look primarily at these three columns.

Now for a graphical representation of the cluster as it will be seen by ScopeSim:

```{code-cell} ipython3
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 8))
plt.plot(src.fields[0]["x"], src.fields[0]["y"], '.')
plt.xlabel("X [arcsec]")
plt.ylabel("Y [arcsec]")
```

## Star Grid and Field

These are two functions that are good to test simulations quickly

```{code-cell} ipython3
from scopesim_templates.stellar import star_field, star_grid

field = star_field(n=400, mmin=15, mmax=25, width=20, height=20, filter_name="Ks")
grid = star_grid(n=400, mmin=15, mmax=25, separation=1 , filter_name="Ks")

plt.figure(figsize=(14, 7))
plt.subplot(121)

size =  np.log10(field.fields[0]["weight"])**2
plt.scatter(field.fields[0]["x"], field.fields[0]["y"], s=size, marker="o")

plt.subplot(122)

size =  np.log10(grid.fields[0]["weight"])**2
plt.scatter(grid.fields[0]["x"], grid.fields[0]["y"], s=size, marker="o")
```

In both cases we generated 400 sources between magnitudes 15 (`mmin`) and 25 (`mmax`).

`star_field` places the stars at random, whereas `star_grid` place them in a regular partern controled by `separation` distance.

The size of the simbols illustrate the magnitudes of the stars

## Stars

The core of the `stellar` module is however the `star` function which can create any field according to the user needs

In this case we generate a stellar field following a 2D gaussian distribution with a star of every type in the pickles stellar library

```{code-cell} ipython3
from scopesim_templates.stellar import stars
from spextra.database import SpecLibrary

lib = SpecLibrary("pickles")

spectypes = lib.template_names
nstars = len(spectypes)

x = np.random.randn(nstars) * 10
y = np.random.randn(nstars) * 10
mags = np.linspace(10, 20, nstars)

src = stars(filter_name="J", amplitudes=mags, x=x, y=y, spec_types=spectypes, library="pickles")

src.fields
```
