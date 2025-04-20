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


class ClusterUtilSingletons:
    """Singletons that usually require the internet."""

    _mamajek = None
    _pickles_ms_v = None
    _mamajek_pickles = None
    _f_mass2mv = None
    _f_mass2idx = None

    def fill_singletons(self):
        if self._mamajek is not None:
            return

        self._mamajek = Table.read(DIRNAME / "mamajek_alt.dat", format="ascii.fixed_width")

        # Catalog 2 contains only main sequence anyway
        tbl = pyckles.SpectralLibrary("pickles").catalog[2].data
        self._pickles_ms_v = tbl["name"][tbl["metalicity"] == "normal"]

        # Only include Spectral Types for which Pickles has an entry
        self._mamajek_pickles = self._mamajek[[spt in self._pickles_ms_v for spt in self._mamajek["SpT"]]]

        self._f_mass2mv = interp1d(
            self._mamajek_pickles["Msun"],
            self._mamajek_pickles["Mv"],
            kind=1,
            # bounds_error=False,
            fill_value="extrapolate",
        )

        self._f_mass2idx = interp1d(
            self._mamajek_pickles["Msun"],
            range(len(self._mamajek_pickles)),
            kind=0,
            bounds_error=False,
            fill_value=(len(self._mamajek_pickles) - 1, 0),
        )

    @property
    def mamajek(self):
        self.fill_singletons()
        return self._mamajek

    @property
    def pickles_ms_v(self):
        self.fill_singletons()
        return self._pickles_ms_v

    @property
    def mamajek_pickles(self) -> np.ndarray:
        self.fill_singletons()
        return self._mamajek_pickles

    @property
    def f_mass2mv(self):
        self.fill_singletons()
        return self._f_mass2mv

    @property
    def f_mass2idx(self):
        self.fill_singletons()
        return self._f_mass2idx


CLUSTER_UTIL_SINGLETONS = ClusterUtilSingletons()


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

    idx = CLUSTER_UTIL_SINGLETONS.f_mass2idx(mass).astype(int)
    return CLUSTER_UTIL_SINGLETONS.mamajek_pickles["SpT"][idx]


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
        return CLUSTER_UTIL_SINGLETONS.f_mass2mv(np.asarray(mass)).round(3)
    # Round to match interpolation precision, float to not have array.
    return float(CLUSTER_UTIL_SINGLETONS.f_mass2mv(mass).round(3))


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

    lum_types = CLUSTER_UTIL_SINGLETONS.pickles_ms_v[CLUSTER_UTIL_SINGLETONS.pickles_ms_v.startswith(spt[0])]
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
