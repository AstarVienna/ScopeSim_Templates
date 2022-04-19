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


def function_call_str(func, args):
    func_str = ".".join([func.__module__, func.__name__])
    args_str = ", ".join(["{}={}".format(key, args[key]) for key in args])
    func_call = "{}({})".format(func_str, args_str)

    return func_call


def make_img_wcs_header(ra, dec, pixel_scale, image_size):
    """
    Create a WCS header for an image

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
    wcs.wcs.crpix = [x // 2, y // 2]
    wcs.wcs.cdelt = np.array([-pixel_scale / 3600, pixel_scale / 3600])
    wcs.wcs.crval = [coords.ra.value, coords.dec.value]

    wcs.wcs.cunit = [u.deg, u.deg]

    return wcs.to_header()


def make_cube_wcs_header():
    """
    TODO: Think if we need it
    """
    return NotImplementedError
