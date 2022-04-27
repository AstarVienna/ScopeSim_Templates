import warnings

import numpy as np
from astropy import units as u
from astropy.io import fits

from spextra import Spextrum

from scopesim_templates.rc import Source, __config__
from scopesim_templates.utils.general_utils import function_call_str, make_img_wcs_header
from ..misc.misc import uniform_source


__all__ = ["lamp",
           "flat_field",
           "empty_sky"]


def lamp(waves, fwhm, fluxes):
    """
    Simple lamp function that creates an homogenous source with a spectra of emission lines (superposed to a very faint
    continuum)

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
    params = locals()
    print("make sure to turn off the appropriate ScopeSim effects during the simulation")
    params["object"] = "lamp"
    params["function_call"] = function_call_str(lamp, params)

    if isinstance(waves, u.Quantity) is False:
        waves = waves*u.AA
    if isinstance(fwhm, u.Quantity) is False:
        fwhm = fwhm*u.AA

    waves = waves.to(u.AA)
    w_min, w_max = np.min(waves), np.max(waves)
    sp = Spextrum.flat_spectrum(amplitude=40, waves=[w_min.value, w_max.value]*u.AA)  # A very faint spextrum to add the lines
    sp = sp.add_emi_lines(center=waves, fwhm=fwhm, flux=fluxes)
    src = uniform_source(sed=sp)
    src.meta.update(params)

    return src


def flat_field(temperature=5000, amplitude=0*u.ABmag, filter_curve="V", extend=60):
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
        any filter curve available for ``spextra``, used for scaling


    Returns
    -------
    src: Source

    """
    params = locals()
    print("make sure to turn off the appropriate ScopeSim effects during the simulation")
    params["object"] = "flat_field"
    params["function_call"] = function_call_str(flat_field, params)

    sp = Spextrum.black_body_spectrum(temperature=temperature, amplitude=amplitude, filter_curve=filter_curve)
    src = uniform_source(sed=sp, amplitude=amplitude, filter_curve=filter_curve, extend=extend)
    src.meta.update(params)

    return src


def empty_sky():
    """
    Returns an empty source so that instrumental fluxes can be simulated

    Returns
    -------
    sky : Source

    """
    params = {"function_call": function_call_str(empty_sky, {}),
              "object": "empty sky"}

    sky = Source(lam=np.array([__config__["!spectral.wave_min"],
                               __config__["!spectral.wave_max"]]),
                 spectra=np.array([0, 0]), x=[0], y=[0], ref=[0], weight=[0])
    sky.meta.update(params)

    return sky
