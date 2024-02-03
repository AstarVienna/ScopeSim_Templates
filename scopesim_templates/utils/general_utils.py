# -*- coding: utf-8 -*-
"""TBA."""

import logging
import functools
from inspect import signature

import numpy as np

from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from synphot import SourceSpectrum, Empirical1D, ConstFlux1D

# TODO: Change this when scopesim is ready RA0, DEC0 = -10, 10
RA0, DEC0 = 0, 0


def hdu_to_synphot(hdu):

    wave = hdu.data["wavelength"]
    wave_unit = u.Unit(hdu.header["TUNIT1"])
    flux = hdu.data["flux"]
    flux_unit = u.Unit(hdu.header["TUNIT2"])

    assert len(wave) > 3, "At least 4 wavelength points are required."
    spec = SourceSpectrum(Empirical1D, points=wave*wave_unit,
                          lookup_table=flux*flux_unit)

    return spec


def vega_spectrum(mag=0):
    vega = SourceSpectrum.from_vega(cache=True)
    return vega * 10**(-0.4 * mag)


def st_spectrum(mag=0):
    return SourceSpectrum(ConstFlux1D, amplitude=mag*u.STmag)


def ab_spectrum(mag=0):
    return SourceSpectrum(ConstFlux1D, amplitude=mag*u.ABmag)


def _get_kwargs_and_defaults(args, kwargs, defaults):
    # BUG: This sometimes raises StopIteration, find out why and catch it!
    default_keys = iter(defaults)
    for arg in args:
        yield next(default_keys), arg
    for key in default_keys:
        yield key, kwargs.get(key, defaults[key].default)


def function_call_str(func, args=None, kwargs=None) -> str:
    """Generate function call signature string including default values."""
    args = args or ()
    kwargs = kwargs or {}
    defaults = signature(func).parameters
    call_kwargs = dict(_get_kwargs_and_defaults(args, kwargs, defaults))
    kwargs_repr = ", ".join(f"{k}={v!r}" for k, v in call_kwargs.items())
    func_str = f"{func.__module__}.{func.__qualname__}"
    call_str = f"{func_str}({kwargs_repr})"
    return call_str


def add_function_call_str(func):
    """Add constructor signature to meta of returned source. Decorator."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        call_str = function_call_str(func, args, kwargs)
        logging.debug("Calling %s ...", call_str)
        src = func(*args, **kwargs)
        src.meta["object"] = func.__name__
        src.meta["function_call"] = call_str
        return src
    return wrapper


def make_img_wcs_header(ra, dec, pixel_scale, image_size):
    """
    Create a WCS header for an image.

    ra : str, float
        "hh:mm:ss.s" or deg for the center of the image
    dec : str, float
        "dd:mm:ss.s or deg for the center of the image
    pixel_scale : float
        arcsecs
    image_size : tuple
        x, y where x, y are integers

    """
    if isinstance(ra, str):  # just assume is in the classical formats
        ra_unit = u.hourangle
    else:
        ra_unit = u.deg

    if isinstance(pixel_scale, u.Quantity):
        pixel_scale = pixel_scale.to(u.arcsec).value

    coords = SkyCoord(ra, dec, unit=(ra_unit, u.deg))
    x, y = image_size

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = [u.deg, u.deg]
    wcs.wcs.crpix = [(x + 1) / 2, (y + 1) / 2]
    wcs.wcs.cdelt = np.array([-pixel_scale / 3600, pixel_scale / 3600])
    wcs.wcs.crval = [coords.ra.value, coords.dec.value]

    wcs.wcs.cunit = [u.deg, u.deg]

    return wcs.to_header()


def make_cube_wcs_header():
    """TODO: Think if we need it."""
    raise NotImplementedError()
