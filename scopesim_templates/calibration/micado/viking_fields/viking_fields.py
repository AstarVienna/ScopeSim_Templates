import os
from os import path as pth
from copy import deepcopy

import numpy as np
from astropy.io import ascii
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import pyckles

from ....stellar import stars
from ....extragalactic import elliptical
from ....rc import ter_curve_utils as tcu


DATA_DIR = pth.abspath(pth.dirname(__file__))

def load_star_sources(fname="cat_illum.dat"):
    tbl = ascii.read(pth.join(DATA_DIR, fname))
    ra0, dec0, mags = np.average(tbl["RA"]), np.average(tbl["Dec"]), tbl["Hmag"]
    x, y = 3600 * (tbl["RA"] - ra0), 3600 * (tbl["Dec"] - dec0)

    src = stars(filter_name="H", flux=mags, spec_types="A0V", x=x, y=y)

    return src


def load_galaxy_sources(fname="faintgal_1.dat", pixel_scale=0.0015):
    """
    Parameters
    ----------
    fname : str
    pixel_scale : float
        [arcsec]

    Returns
    -------
    srcs : list of Source

    """
    tbl = ascii.read(pth.join(DATA_DIR, fname))

    H_mag = tbl["H_mag"].data.astype(float)         # [Vega mag]
    Reff = tbl["Reff"].data.astype(float)           # [mas]
    gal_type = tbl["type"].data.astype(float)       # []

    # Or should there be a range > involved? E.g. Spiral 1-2, Ellip 3-5?
    # Ellip type 0, spiral type 1
    # Do spirals first
    n_gals = len(gal_type)
    ns = np.random.random(n_gals)   # Sersic index
    ns[gal_type==1] += 1            # Spirals n=[1..2]
    ns[gal_type==0] *= 2            # Ellips n=[3..5] => ([0..1] * 2) + 3
    ns[gal_type==0] += 3

    scale_factors = 10**(-0.4*H_mag)
    r_effs = 1e-3 * Reff                            # [arcsec]
    xs, ys = np.random.random((2, n_gals)) * 60     # [arcsec]
    angles = np.random.random(n_gals) * 360         # [deg]
    ellipticitys = np.ones(n_gals)
    e_max = 0.85                        # 0 = circular
    ellipticitys = e_max * np.random.random(n_gals)
    ellipticitys[gal_type == 0] = 0
    widths = 2 * r_effs

    brown_lib = pyckles.SpectralLibrary("brown", return_style="synphot")
    spectrum = brown_lib["NGC_0584"]
    spectrum = tcu.scale_spectrum(spectrum, filter_name="H", amplitude=0*u.mag)
    spectra = [sf * deepcopy(spectrum) for sf in scale_factors]

    srcs = []
    for i in range(n_gals):
        src = elliptical(r_eff=r_effs[i],
                         pixel_scale=pixel_scale,
                         filter_name="H",
                         amplitude=1,
                         spectrum=spectra[i],
                         rescale_spectrum=False,
                         n=ns[i],
                         angle=angles[i],
                         ellipticity=ellipticitys[i],
                         width=widths[i],
                         height=widths[i],
                         normalization="total",
                         )
        src.shift(dx=xs[i], dy=ys[i])
        print(i, r_effs[i], widths[i])
        srcs += [src]

    return srcs


