# -*- coding: utf-8 -*-
"""WCU laser spectrum for METIS."""

import numpy as np
from astropy import units as u
from astropy.io import fits

from spextra import Spextrum

from ..utils.general_utils import add_function_call_str
from ..rc import Source, im_plane_utils


@add_function_call_str
def laser_spectrum(
        centers,
        fwhms,
        fluxes,
        img_size=(12, 12),
        waves=None,
        specdict=None,
):
    """
    Create laser calibration spectrum source for METIS.

    Parameters
    ----------
    centers : TYPE
        Emission line centers.
    fwhms : TYPE
        Emission line FWHMs.
    fluxes : TYPE
        Emission line total fluxes.
    img_size : TYPE, optional
        FOV size in arcsec, plus a little extra. The default is (12, 12).
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
    spec = Spextrum.emission_line_spectrum(centers.round(5), fwhms, fluxes,
                                           waves=waves)
    src = Source(image_hdu=field, spectra=spec)
    return src


def laser_spectrum_lm(specdict, n_tunable: int = 20, **kwargs):
    """
    Calibration laser for METIS LM bands.

    Create a fixed red and blue laser as well as a variable number of tunable
    laser lines in between.

    Parameters
    ----------
    specdict : dict-like
        Dict-like structure containing the keys "wave_min", "wave_max",
        "spectral_bin_width" and "wave_unit", which is used to construct the
        wavelength axis. Usually this means ``cmd["!SIM.spectral"]``.
    n_tunable : int, optional
        Number of emission lines created for the tunable laser between 4.68 and
        4.78 um. Passing 1 will create a single line at 4.73 um.
        The default is 20.
    **kwargs : TYPE
        Any other keyword arguments passed to ``laser_spectrum``.
        In particular, `centers`, `fwhms` and `fluxes` will override the
        default values and any lines created by `n_tunable`. If those are
        overridding, make sure the shapes of those three arguments match.

    Returns
    -------
    src : Source

    Notes
    -----
    The default values for the emission lines are 3.39 and 5.263 um for the
    fixed lasers, with FWHM of 8E-8 and 6E-5 um respectively. The various lines
    created by the tunable laser all use the same FWHM of 8E-7 um. The fluxes
    of all lines are scaled to create lines of roughly equal height.

    """
    if n_tunable == 1:
        tunable_cen = [4.73] * u.um
    else:
        tunable_cen = np.linspace(4.68, 4.78, n_tunable) * u.um

    defaults = {
        "centers": [3.39, *tunable_cen.value, 5.263] * u.um,
        "fwhms": [8E-8, *[8E-7]*len(tunable_cen), 6E-5] * u.um,
        "fluxes": [1, *[7]*len(tunable_cen), 500],
        "specdict": specdict
    } | kwargs
    return laser_spectrum(**defaults)


def laser_spectrum_n(specdict, **kwargs):
    """
    Calibration laser for METIS N band.

    Create equal emission lines at 9, 10, 11 and 12 um.

    Parameters
    ----------
    specdict : dict-like
        Dict-like structure containing the keys "wave_min", "wave_max",
        "spectral_bin_width" and "wave_unit", which is used to construct the
        wavelength axis. Usually this means ``cmd["!SIM.spectral"]``.
    **kwargs : TYPE
        Any other keyword arguments passed to ``laser_spectrum``.

    Returns
    -------
    src : Source

    Notes
    -----
    The default of 1E-8 um is used as a best-guess for the line's FWHMs.
    The relative fluxes are scaled to produce lines of roughly equal height.
    """
    defaults = {
        "centers": [9, 10, 11, 12] * u.um,
        "fwhms": 4 * [1E-8] * u.um,
        "fluxes": [1.3, 1.2, 1.1, 1],
        "specdict": specdict
    } | kwargs
    return laser_spectrum(**defaults)
