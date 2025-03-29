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

# Extragalactic

A few functions that can also helpt to create extragalactic objects

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

from scopesim_templates.extragalactic import galaxy
```

The function `galaxy` will generate a sersic profile with user provided parameters. The function can accept any sed in the [speXtra](https://spextra.readthedocs.io/en/latest/) package.

It must be noted that the pixel scale is not related to the pixel scale of the final simulation but to the pixel scale (in arcsec) of the generated image. It is recommended to be fine enough to well sample the profile of the galaxy.

The SED and the selected filter must overlap but an error will be thrown if they don't.

```{code-cell} ipython3
src = galaxy("kc96/s0", z=0.1, amplitude=17, filter_curve="g", pixel_scale=0.05, r_eff=2.5, n=4, ellip=0.5, theta=45, extend=3)

plt.imshow(np.log10(src.fields[0].data), origin="lower")
```

We can also visualize the spectrum

```{code-cell} ipython3
sp = src.spectra[0]
sp.plot()
```

Finally there is much that can be done with the spectrum itself, like for example, obtaining the magnitude in a different filter

```{code-cell} ipython3
sp.get_magnitude(filter_curve="u", system_name="AB")
```

More can be found in the [speXtra](https://spextra.readthedocs.io/en/latest/) package.
