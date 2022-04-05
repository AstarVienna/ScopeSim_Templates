import numpy as np
from astropy import units as u
from astropy.io import fits
from scopesim_templates.rc import Source, __config__
from scopesim_templates.utils.general_utils import function_call_str, make_img_wcs_header
from spextra import Spextrum


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
        any filter curve available for ``spextra``


    Returns
    -------
    src: Source

    """

    sp = Spextrum.black_body_spectrum(temperature=temperature, amplitude=amplitude, filter_curve=filter_curve)

    if isinstance(amplitude, u.Quantity) is False:
        amplitude = amplitude * u.ABmag

    data = np.ones(shape=(extend, extend))
    header = make_img_wcs_header(ra=0, dec=0, pixel_scale=1, image_size=data.shape)
    hdu = fits.ImageHDU(header=header, data=data)

    src = Source(spectra=sp, image_hdu=hdu)

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
