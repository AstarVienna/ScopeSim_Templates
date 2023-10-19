"""Flat lamp for MICADO."""

import numpy as np
from astropy.io import fits
import astropy.units as u

from ..rc import Source


def flatlamp(
        width=16.384,
        height=16.384,
        pixel_scale=1.024,
        amplitude=1000.0,
        fraction=0.9,
    ):
    """
    A flatlamp.

    Default parameters selected such that:
    - image is 16x16
    - covering a single 4k by 4k detector.

    Parameters
    ----------
    width : float
        [arcsec] width of the 'screen'
    height : float
        [arcsec] width of the 'screen'
    pixel_scale : float
        [arcsec/pix] pixel scale
    amplitude : float
        [None] arbitrary scale factor
    fraction : float
        [None] fraction of flux at the edge

    Returns
    -------
    src: scopesim.Source
    """

    # 16x16 by default
    sx, sy = int(width / pixel_scale), int(height / pixel_scale)
    ax = np.arange(sx)
    ay = np.arange(sy)
    xx, yy = np.meshgrid(ax, ay, sparse=True)

    crpx = (width - 1) / 2
    crpy = (height - 1) / 2

    # * 1. is necessary for the normalization later
    image = - ((xx - crpx) ** 2 + (yy - crpy) ** 2) * 1.

    imin = image.min()
    imax = image.max()
    image = (1. - fraction) * (image - imin) / (imax - imin) + fraction
    image *= amplitude

    hdu = fits.ImageHDU(data=image)

    hdu.header["CRPIX1"] = crpx + 1
    hdu.header["CRPIX2"] = crpy + 1
    hdu.header["CRVAL1"] = 0
    hdu.header["CRVAL2"] = 0
    hdu.header["CDELT1"] = pixel_scale / 3600
    hdu.header["CDELT2"] = pixel_scale / 3600
    hdu.header["CUNIT1"] = "deg"
    hdu.header["CUNIT2"] = "deg"
    hdu.header["CTYPE1"] = "RA---TAN"
    hdu.header["CTYPE2"] = "DEC--TAN"

    # number_of_points should be at least 200 or so
    number_of_points = 230
    sm = 1.2419516582773063
    spectra = sm * np.ones(number_of_points)

    sa = 1.0200000e+02
    se = 3.3502039e+05
    lam = np.logspace(np.log(sa), np.log(se), base=np.e,
                      num=number_of_points) * u.Angstrom

    return Source(
            image_hdu=hdu,
            spectra=spectra,
            lam=lam,
        )
