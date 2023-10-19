import numpy as np
from astropy import units as u
from synphot import BlackBody1D, SourceSpectrum, Empirical1D

from ..rc import Source


def pinhole_mask(x, y, waves, temperature=1500*u.K, sum_factor=1):
    """
    Creates a Source object with point sources at coordinates (x, y) in the slit

    Parameters
    ----------
    x, y : array-like
        [arcsec] Coords relative to centre of FOV (i.e. optical axis)
    waves : array-like, u.Quantity
        [um] Array of wavelength values.
    temperature : float, u.Quantity, optional
        [deg K] Temperature of WCU blackbody spectrum.
    sum_factor : float, optional
        Rescale the spectrum to sum to this value


    Returns
    -------
    src : scopesim.Source


    Examples
    --------
    Pin-holes every 0.5 arcsec along the MICADO long-slit
    ::
        x = np.arange(-1.5, 13.51, 0.5)
        y = np.zeros_like(x)
        waves = np.arange(0.7, 2.5, 0.001) * u.um

        src = pinhole_mask(x, y, waves, sum_factor=9001)

    A grid of pin-holes every 5 arcsec for the MICADO 4mas IMG mode
    ::
        dr = np.arange(-25, 26, 5)      # [arcsec]
        x, y = np.meshgrid(dr, dr)
        x, y = x.flatten(), y.flatten()
        waves = np.arange(0.7, 2.5, 0.001) * u.um

        src = pinhole_mask(x, y, waves, sum_factor=9001)

    """
    if not isinstance(waves, u.Quantity):
        waves *= u.um
    if not isinstance(temperature, u.Quantity):
        temperature *= u.K

    blackbody_spectrum = BlackBody1D(temperature=temperature)
    flux = blackbody_spectrum(waves.to(u.AA))
    flux *= sum_factor / np.trapz(flux)
    spec = SourceSpectrum(Empirical1D, points=waves, lookup_table=flux)
    src = Source(x=x, y=y, ref=np.ones_like(x), weight=np.ones_like(x),
                 spectra=[spec])

    return src
