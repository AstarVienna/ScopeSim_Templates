# -*- coding: utf-8 -*-
"""Contains aux functions for clusters.py."""

from pathlib import Path
from collections.abc import Iterable

import numpy as np
from astropy.table import Table
# from astropy.modeling.functional_models import KingProjectedAnalytic1D
from scipy.interpolate import interp1d

import pyckles


DIRNAME = Path(__file__).parent
MAMAJEK = Table.read(DIRNAME / "mamajek_alt.dat", format="ascii.fixed_width")

# Catalog 2 contains only main sequence anyway
tbl = pyckles.SpectralLibrary("pickles").catalog[2].data
PICKLES_MS_V = tbl["name"][tbl["metalicity"] == "normal"]

# Only include Spectral Types for which Pickles has an entry
MAMAJEK_PICKLES = MAMAJEK[[spt in PICKLES_MS_V for spt in MAMAJEK["SpT"]]]
F_MASS2MV = interp1d(
    MAMAJEK_PICKLES["Msun"],
    MAMAJEK_PICKLES["Mv"],
    kind=1,
    # bounds_error=False,
    fill_value="extrapolate",
)
F_MASS2IDX = interp1d(
    MAMAJEK_PICKLES["Msun"],
    range(len(MAMAJEK_PICKLES)),
    kind=0,
    bounds_error=False,
    fill_value=(len(MAMAJEK_PICKLES) - 1, 0),
)


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
        mass = np.asarray(mass)

    idx = F_MASS2IDX(mass).astype(int)
    return MAMAJEK_PICKLES["SpT"][idx]


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
        return F_MASS2MV(np.asarray(mass)).round(3)
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
