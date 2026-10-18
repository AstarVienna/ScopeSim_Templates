# -*- coding: utf-8 -*-
"""TBA."""

import logging
import functools
from inspect import signature, Parameter

import numpy as np

from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from synphot import SourceSpectrum, Empirical1D, ConstFlux1D

# TODO: Change this when scopesim is ready RA0, DEC0 = -10, 10
RA0, DEC0 = 0, 0


def get_rng(seed=None) -> np.random.Generator:
    """
    Return a NumPy random generator for the given seed.

    All stochastic templates draw their random numbers from the generator
    returned here and never touch NumPy's global random state.

    Parameters
    ----------
    seed : None, int, array_like[int], np.random.SeedSequence or np.random.Generator
        Anything accepted by ``np.random.default_rng``. A ``Generator`` is
        returned unchanged, so a caller can share one stream between several
        calls. If None, ``rc.__config__["!random.seed"]`` is used, which is
        itself None by default (i.e. a fresh realisation on every call).

    Returns
    -------
    rng : np.random.Generator

    """
    if seed is None:
        from ..rc import __config__
        seed = __config__["!random.seed"]
    return np.random.default_rng(seed)


def hdu_to_synphot(hdu):

    wave = hdu.data["wavelength"]
    wave_unit = u.Unit(hdu.header["TUNIT1"])
    flux = hdu.data["flux"]
    flux_unit = u.Unit(hdu.header["TUNIT2"])

    if len(wave) < 4:
        raise ValueError("At least 4 wavelength points are required.")
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


def _get_kwargs_and_defaults(args, kwargs, func):
    """Yield (name, value) for every argument of a call, defaults included."""
    sig = signature(func)
    try:
        bound = sig.bind_partial(*args, **kwargs)
    except TypeError:
        # The call itself will fail with a proper error message; just
        # report what was passed.
        yield from ((f"arg{i}", arg) for i, arg in enumerate(args))
        yield from kwargs.items()
        return

    bound.apply_defaults()
    for name, value in bound.arguments.items():
        kind = sig.parameters[name].kind
        if kind is Parameter.VAR_KEYWORD:
            yield from value.items()
        elif kind is Parameter.VAR_POSITIONAL:
            yield from ((f"{name}[{i}]", arg) for i, arg in enumerate(value))
        else:
            yield name, value


def function_call_str(func, args=None, kwargs=None) -> str:
    """
    Generate function call signature string including default values.

    This is only used to document how a Source was made, so it never raises:
    arguments that don't fit the signature are listed as given.
    """
    args = args or ()
    kwargs = kwargs or {}
    call_kwargs = dict(_get_kwargs_and_defaults(args, kwargs, func))
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


# TODO: This is WET with create_wcs_from_points, try to unifiy those.
#       Also, the whole ra_unit stuff can probably be handled by astropy...
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
        (nx, ny), i.e. the number of pixels along x (NAXIS1) and y (NAXIS2).
        Note this is the *reverse* of a numpy ``array.shape``.

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
