import numpy as np
from astropy.io import ascii
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scopesim_templates.stellar import stars


def load_star_cat(fname="cat_illum.dat"):
    tbl = ascii.read(fname)
    ra0, dec0, mags = np.average(tbl["RA"]), np.average(tbl["Dec"]), tbl["Hmag"]
    x, y = 3600 * (tbl["RA"] - ra0), 3600 * (tbl["Dec"] - dec0)

    src = stars(filter_name="H", flux=mags, spec_types="A0V", x=x, y=y)

    return src


def load_galaxy_cat(fname="faintgal_1.dat"):
    tbl = ascii.read(fname)


