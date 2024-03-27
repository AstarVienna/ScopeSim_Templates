# -*- coding: utf-8 -*-
"""WCU laser spectrum for METIS."""

import numpy as np
from astropy import units as u
from astropy.io import fits

from spextra import Spextrum

from ..rc import Source, im_plane_utils


def laser_spectrum(img_size=(12, 12),
                   centers=[3.39, 4.73, 5.263]*u.um,
                   fwhms=[8E-8, 8E-7, 6E-5]*u.um,
                   fluxes=[1, 7, 500],
                   waves=None,
                   specdict=None,
):
    """
    Create laser calibration spectrum source for METIS.

    Parameters
    ----------
    img_size : TYPE, optional
        FOV size in arcsec, plus a little extra. The default is (12, 12).
    centers : TYPE, optional
        Emission line centers. The default is [3.39, 4.73, 5.263]*u.um.
    fwhms : TYPE, optional
        Emission line FWHMs. The default is [8E-8, 8E-7, 6E-5]*u.um.
    fluxes : TYPE, optional
        Emission line total fluxes. The default is [1, 7, 500] to create lines
        of roughly equal amplitude.
    waves : array-like, optional
        Input wavelength axis. Either `waves` or `specdict` must be given. If
        both are given, `waves` takes precedence.The default is None.
    specdict : Mapping, optional
        Dict-like structure containing the keys "wave_min", "wave_max",
        "spectral_bin_width" and "wave_unit", which is used to construct the
        wavelength axis. Either `waves` or `specdict` must be given. If both
        are given, `waves` takes precedence. The default is None.

    Returns
    -------
    src : Source

    """
    if waves is None:
        if specdict is None:
            raise ValueError("Either waves or specdict must be given.")
        waves = np.arange(
            specdict["wave_min"],
            specdict["wave_max"],
            specdict["spectral_bin_width"]
        ) * u.Unit(specdict["wave_unit"])

    # Taken from micado.spectral_calibration.line_list:
    dw = 0.5 * img_size[0] / 3600         # input needed in [deg]
    dh = 0.5 * img_size[1] / 3600         # input needed in [deg]
    pixel_scale = 0.1 * dw

    hdr = im_plane_utils.header_from_list_of_xy(
        [-dw, dw], [-dh, dh], pixel_scale)  # [deg]
    hdr["SPEC_REF"] = 0
    im = np.ones((hdr["NAXIS1"], hdr["NAXIS2"]))
    field = fits.ImageHDU(header=hdr, data=im)
    spec = Spextrum.emission_line_spectrum(centers, fwhms, fluxes, waves=waves)
    src = Source(image_hdu=field, spectra=spec)
    return src
