from pathlib import Path

import numpy as np
from scipy import signal
from astropy.io import ascii as ioascii
from astropy.io import fits
from astropy import units as u
from synphot import SourceSpectrum, Empirical1D
from synphot.units import PHOTLAM

from ..rc import Source, im_plane_utils

DATA_DIR = Path(__file__).parent / "data"


def line_list(unit_flux=1*PHOTLAM,
              dwave=0.001,
              smoothing_fwhm=0.003,
              width=15.,
              height=1.,
              filename="masterlinelist_2.1.txt",
              ):
    """
    Make a ScopeSim.Source object for the calibration lamp line list.

    [PHOTLAM] = [ph s-1 cm-2 AA-1]

    Parameters
    ----------
    unit_flux : astropy.Quantity
        Default 1*PHOTLAM.
        Flux conversion factor for a relative intensity value of 1.
        The unit can be overridden with any equivalent astropy Quantity.
    dwave : float
        Bin width in units of input file
    smoothing_fwhm : int, float
        FWHM of smoothing kernel in units of ``dwave`` (i.e. input file)
    filename : str
        Name of line list file found in ``micado/data``

    Returns
    -------
    line_list_src : scopesim.Source

    Examples
    --------
    >>> from scopesim_templates.micado import spectral_calibrations as mic_spec
    >>> from synphot.units import PHOTLAM
    >>> from astropy import units as u
    >>>
    >>> # Use all defaults
    >>> src = mic_spec.line_list()
    >>>
    >>> # Make line list brighter
    >>> src = mic_spec.line_list(unit_flux=2000)        # [PHOTLAM] by default
    >>> src = mic_spec.line_list(unit_flux=5*PHOTLAM)
    >>> src = mic_spec.line_list(unit_flux=300*u.Unit("ph s-1 m-2 um-1"))
    >>>
    >>> # Smooth with a gaussian kernel where FWHM is 3 bins wide
    >>> src = mic_spec.line_list(smoothing_fwhm=3)

    """
    if not isinstance(unit_flux, u.Quantity):
        unit_flux *= PHOTLAM

    # calculate smoothing kernel size
    if do_smoothing := smoothing_fwhm is not None:
        sigma = smoothing_fwhm / dwave / 2.35
        ksize = int(8 * sigma + 1)
    else:
        ksize = 2

    # import line list and pad with zeros
    # (including smoothing kernel space at the edges)
    wave, flux = import_line_spectrum(filename, dwave, pad=ksize)

    if do_smoothing:
        kernel = signal.windows.gaussian(M=ksize, std=sigma)
        kernel /= kernel.sum()
        flux = signal.convolve(flux, kernel, mode="same")

    # make spectrum
    assert len(wave) > 3, "At least 4 wavelength points are required."
    spec = SourceSpectrum(Empirical1D, points=wave * u.nm,
                          lookup_table=flux * unit_flux)

    # make unity image that covers MICADO FOV of +/- 30" arcsec
    dw = 0.5 * width / 3600         # input needed in [deg]
    dh = 0.5 * height / 3600         # input needed in [deg]
    pixel_scale = 0.1 * dw

    hdr = im_plane_utils.header_from_list_of_xy([-dw, dw], [-dh, dh],
                                                pixel_scale)  # [deg]
    hdr["SPEC_REF"] = 0
    im = np.ones((hdr["NAXIS1"], hdr["NAXIS2"]))
    field = fits.ImageHDU(header=hdr, data=im)

    line_list_src = Source(image_hdu=field, spectra=[spec])

    return line_list_src


def import_line_spectrum(filename, dwave=0.0001, pad=10):
    """
    Read in and pad the line list into a spectrum.

    Parameters
    ----------
    filename : str
        Name of line list file found in ``micado/data``

    dwave : float
        Bin width in units of input file

    pad : int
        Padding in wavelength range to add to allow enough space for smoothing
        function kernel

    Returns
    -------
    wave_filled : np.ndarray
        wavelength vector in units of whatever is in the input file
    flux_filled : np.ndarray
        flux vector in units of whatever is in the input file

    """
    line_tbl = ioascii.read(DATA_DIR / filename)

    flux = line_tbl["Relative_Intensity"].data
    wave = line_tbl["Wavelength"].data
    wave = np.round(wave / dwave) * dwave       # round to level of dwave
    padwave = pad * dwave
    wave_filled = np.arange(wave.min() - padwave, wave.max() + padwave, dwave)
    flux_filled = np.zeros_like(wave_filled)

    indices = np.searchsorted(wave_filled, wave)
    np.add.at(flux_filled, indices, flux)

    return wave_filled, flux_filled
