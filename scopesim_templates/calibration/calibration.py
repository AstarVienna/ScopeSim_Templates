# -*- coding: utf-8 -*-
"""Contains simple source functions for testing and calibration.

Simple templates that could be used to simulate calibration frames.
Make sure to turn off the corresponding effects during the simulation.
"""

import warnings

import numpy as np
from astropy import units as u

from spextra import Spextrum, Passband

from ..rc import Source, __config__
from ..utils.general_utils import add_function_call_str
from ..misc.misc import uniform_source

__all__ = [
    "lamp",
    "flat_field",
    "empty_sky",
]


@add_function_call_str
def lamp(waves, fwhm, fluxes) -> Source:
    """
    Create a homogenous source with a spectrum of emission lines.

    Emission line spectrum is superposed to a very faint continuum.

    Parameters
    ----------
    waves : np.array, u.Quantity
        a list or array of wavelengths
    fwhm : np.array, u.Quantity
       a list or array with the fwhm of the lines
    fluxes : np.array, u.Quantity
       a list or array with the fluxes of the lines

    Return
    ------
    src : Source

    Examples
    --------
    Create a ``Source`` lamp::

        >>> import numpy as np
        >>> import astropy.units as u
        >>> from scopesim_templates.calibration import lamp
        >>>
        >>> R = 20000   # Resolution of the instrument
        >>> waves = np.arange(1e4, 2.7e4, 100)
        >>> fwhm = 2.6*(waves/R)*np.ones(shape=waves.shape)  # Nyquist sampled
        >>> flux = 1e-10*np.ones(waves.shape)
        >>>
        >>> src = lamp(waves=waves, fwhm=fwhm, fluxes=flux)

    """
    # I don't understand this warning but it probably has a reason...
    warnings.warn("Make sure to turn off the appropriate ScopeSim effects "
                  "during the simulation", RuntimeWarning)

    waves <<= u.AA
    fwhm <<= u.AA

    w_minmax = [waves.min().value, waves.max().value] * u.AA
    # A very faint spextrum to add the lines
    spec = Spextrum.flat_spectrum(
        amplitude=40 * u.ABmag,
        waves=w_minmax,
    )
    spec = spec.add_emi_lines(center=waves, fwhm=fwhm, flux=fluxes)
    src = uniform_source(sed=spec)
    # src.meta.update({"object": "lamp"})

    return src


@add_function_call_str
def flat_field(
    temperature: u.Quantity[u.K] = 5000 * u.K,
    amplitude: u.Quantity[u.ABmag] = 0 * u.ABmag,
    filter_curve: str | Passband = "V",
    extend=60,
) -> Source:
    """
    Create a flat-field source to be used in ``ScopeSim``.

    The spectral shape is given by a Black Body Spectrum with a temperature set
    by the user, so it can also be used in spectroscopy with enough realism.

    Flats-fields usually also contain an illumination pattern. They might be
    implemented in ScopeSim.effects.
    TODO: Investigate if that belongs eventually here.

    Default values are just a wild guess. We need to find more realistic ones.

    Parameters
    ----------
    temperature: ``astropy.Quantity``
        [Kelvin] Temperature of the lamp
    amplitude: ``astropy.Quantity``
        [u.Quantity] amplitude of the lamp in u.ABmag, u.mag (vega) or u. STmag
    filter_curve: str
        any filter curve available for ``spextra``, used for scaling


    Returns
    -------
    src: Source
    """
    # I don't understand this warning but it probably has a reason...
    warnings.warn("Make sure to turn off the appropriate ScopeSim effects "
                  "during the simulation", RuntimeWarning)

    spec = Spextrum.black_body_spectrum(
        temperature=temperature,
        amplitude=amplitude,
        filter_curve=filter_curve,
    )
    src = uniform_source(
        sed=spec,
        amplitude=amplitude,
        filter_curve=filter_curve,
        extend=extend,
    )
    # src.meta.update({"object": "flat_field"})

    return src


@add_function_call_str
def empty_sky() -> Source:
    """
    Return an empty source so that instrumental fluxes can be simulated.

    Returns
    -------
    sky : Source

    """
    sky = Source(lam=np.array([__config__["!spectral.wave_min"],
                               __config__["!spectral.wave_max"]]),
                 spectra=np.array([0, 0]), x=[0], y=[0], ref=[0], weight=[0])
    # src.meta.update({"object": "empty_sky"})

    return sky
