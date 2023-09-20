# -*- coding: utf-8 -*-
"""Contains aux functions for clusters.py."""

from pathlib import Path
from collections.abc import Iterable

import numpy as np
from astropy.io.ascii import read as read_ascii
# from astropy.modeling.functional_models import KingProjectedAnalytic1D
from scipy.interpolate import interp1d

import pyckles


DIRNAME = Path(__file__).parent
MAMAJEK = read_ascii(DIRNAME / "mamajek_alt.dat", delimiter="|")
F_MASS2MV = interp1d(MAMAJEK["Msun"], MAMAJEK["Mv"], 1)
F_MASS2IDX = interp1d(MAMAJEK["Msun"], range(len(MAMAJEK["Msun"])), 0)
PICKLES = pyckles.SpectralLibrary("pickles")

tbl = PICKLES.catalog[1].data
mask_evol = tbl["evolution"] == 5
mask_metal = tbl["metalicity"] == "normal"
PICKLES_MS_V = tbl["name"][mask_evol * mask_metal]


def mass2spt(mass):
    """
    Find the closest spectral type for a given stellar mass.

    Parameters
    ----------
    mass : float or iterable of floats
        Star mass in units of solar mass.

    Returns
    -------
    spt : str or list of str
        Spectral type or list of spectral types.

    Notes
    -----
    `mass` should not be passed as a quantity.

    """
    if isinstance(mass, Iterable):
        # recursive
        return [mass2spt(m) for m in mass]

    idx = F_MASS2IDX(mass).astype(int)
    return MAMAJEK["SpT"][idx]


def mass2absmag(mass):
    """
    Interpolate the absolute V magnitude (Mv) for a given stellar mass.

    Parameters
    ----------
    mass : float or iterable of floats
        Star mass in units of solar mass.

    Returns
    -------
    absmag : float or list of floats
        Absolute magnitude or magnitudes.

    Notes
    -----
    `mass` should not be passed as a quantity.

    """
    if isinstance(mass, Iterable):
        # recursive
        return [mass2absmag(m) for m in mass]

    # Round to match interpolation precision, float to not have array.
    return float(F_MASS2MV(mass).round(3))


def closest_pickles(spt):
    """
    Retrieve closest spectral type available in the loaded Pyckles catalog.

    Parameters
    ----------
    spt : str or sequence of str
        Spectral type or types.

    Returns
    -------
    closest_pickle : str or list of str
        Closest spectral type or types in the catalog.

    """
    if isinstance(spt, Iterable) and not isinstance(spt, str):
        # recursive
        return [closest_pickles(s) for s in spt]

    to_strip = "OBAFGKMIV "

    lum_types = PICKLES_MS_V[PICKLES_MS_V.startswith(spt[0])]
    assert lum_types.strip(to_strip).isdigit().all()
    # The join is necessary to correctly parse 'M25V' as 2.5 and not 25...
    num_only = np.char.join(".", lum_types.strip(to_strip))
    classes = np.array(num_only).astype(float)

    idx = np.abs(classes - float(spt.strip(to_strip))).argmin()
    closest_pickle = lum_types[idx]

    return closest_pickle


# Unused function, commented out to not annoy pylint...
# def king_distribution(n, r_core, r_tidal):
#     king = KingProjectedAnalytic1D(r_core=r_core, r_tidal=r_tidal)
#     y = king(np.linspace(0, r_tidal, 100))

#     return None


def gaussian_distribution(n, fwhm, seed=None):
    if isinstance(seed, int):
        np.random.seed(seed)
    x, y = np.random.normal(loc=0, scale=fwhm/2.35, size=(2, n))
    return x, y
