from os import path as p

import numpy as np
from scipy import signal
from astropy.io import ascii as ioascii
from astropy.io import fits
from astropy import units as u
from synphot import SourceSpectrum, Empirical1D
from synphot.units import PHOTLAM

from scopesim import Source
from scopesim.optics import image_plane_utils as ipu
from ..rc import ter_curve_utils as tcu

DATA_DIR = p.join(p.dirname(__file__), "data")


def line_list(unit_flux=1*PHOTLAM,
              dwave=0.0001,
              smoothing_fwhm=5,
              filename="masterlinelist_2.1.txt"):
    """
    Make a ScopeSim.Source object for the calibration lamp line list

    [PHOTLAM] = [ph s-1 cm-2 AA-1]

    Parameters
    ----------
    unit_flux : astropy.Quantity
        Default 1*PHOTLAM.
        Flux conversion factor for a relative intensity value of 1.
        The unit can be overridden with any equivalent astropy Quantity.
    dwave : float
        [nm] Bin width
    smoothing_fwhm : int, float
        Default 5 bins. FWHM of smoothing kernel in units of number of bins
        (bin width = dwave)
    filename : str
        Name of line list file found in ``micado/data``

    Returns
    -------
    line_list_src : scopesim.Source

    Examples
    --------
    ::

        from scopesim_templates.micado import spectral_calibrations as mic_spec
        from synphot.units import PHOTLAM
        from astropy import units as u

        # Use all defaults
        src = mic_spec.line_list()

        # Make line list brighter
        src = mic_spec.line_list(unit_flux=2000)        # [PHOTLAM] by default
        src = mic_spec.line_list(unit_flux=5*PHOTLAM)
        src = mic_spec.line_list(unit_flux=300*u.Unit("ph s-1 m-2 um-1))

        # Smooth with a gaussian kernel where FWHM is 3 bins wide
        src = mic_spec.line_list(smoothing_fwhm=3)

    """
    if not isinstance(unit_flux, u.Quantity):
        unit_flux *= PHOTLAM

    # import line list and pad with zeros
    wave, flux = import_line_spectrum(filename)
    if smoothing_fwhm is not None and isinstance(smoothing_fwhm, (int, float)):
        sigma = smoothing_fwhm / 2.35
        kernel = signal.windows.gaussian(M=int(8 * sigma + 1), std=sigma)
        kernel /= kernel.sum()
        flux = signal.convolve(flux, kernel, mode="same")

    # make spectrum
    spec = SourceSpectrum(Empirical1D, points=wave * u.nm,
                          lookup_table=flux * unit_flux)

    # make unity image that covers MICADO FOV of +/- 30" arcsec
    hw = 30./3600         # input needed in [deg]
    hdr = ipu.header_from_list_of_xy([-hw, hw], [-hw, hw], 1./3600)      # [deg]
    hdr["SPEC_REF"] = 0
    im = np.ones((hdr["NAXIS2"], hdr["NAXIS2"]))
    field = fits.ImageHDU(header=hdr, data=im)

    line_list_src = Source(image_hdu=field, spectra=[spec])

    return line_list_src


def import_line_spectrum(filename):
    """
    Read in and pad the line list into a spectrum

    Parameters
    ----------
    filename : str
        Name of line list file found in ``micado/data``

    Returns
    -------
    wave_filled : np.ndarray
        wavelength vector in units of whatever is in the input file
    flux_filled : np.ndarray
        flux vector in units of whatever is in the input file

    """
    line_tbl = ioascii.read(p.join(DATA_DIR, filename))
    dwave = 0.0001
    wave = line_tbl["Wavelength"]
    flux = line_tbl["Relative_Intensity"]
    i_bin = ((wave - wave[0]) / dwave).astype(int)
    wave_filled = np.arange(wave[0], wave[-1] + dwave, dwave)
    flux_filled = np.zeros(i_bin[-1] + 1)
    # .. todo: warning, this only works if dwave is smaller than the minimum
    #          distance between lines
    flux_filled[i_bin] += flux

    return wave_filled, flux_filled
