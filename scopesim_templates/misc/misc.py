import numpy as np

from ..rc import Source, __config__

import astropy.units as u
from spextra import Spextrum


def flat_field(temperature=5000, amplitude=0*u.ABmag, filter_curve="V"):
    """
    This function creates a flat-field source to be used in ScopeSim

    The spectral shape is given by a Black Body Spectrum with a temperature set by the user, so it can
    also be used in spectroscopy with enough realism.

    Flats-fields usually also contain an illumination pattern. They might be implemented in ScopeSim.effects.
    TODO: Investigate if that belongs eventually here.

    Default values are just a wild guess. We need to find more realistic ones.

    Parameters
    ----------
    temperature: float
        [Kelvin] Temperature of the lamp
    amplitude: u.Quantity
        [u.Quantity] amplitude of the lamp in u.ABmag, u.mag (vega) or u. STmag
    filter_curve: str
        any filter curve available for ``spextra``


    Returns
    -------
    src: Source

    """

    sp = Spextrum.black_body_spectrum(temperature=temperature, amplitude=amplitude, filter_curve=filter_curve)

    src = Source(lam=np.array([__config__["!spectral.wave_min"],
                               __config__["!spectral.wave_max"]]),
                 spectra=sp, x=[0], y=[0], ref=[0], weight=[0])
    return src

