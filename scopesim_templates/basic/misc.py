import numpy as np

from ..rc import Source, __config__
from ..utils.general_utils import function_call_str
from astropy.io import fits

from spextra import Spextrum


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


def flat_field(x_size, y_size, counts,  temperature):
    """
    Creates a flat-field image to use in the simulator

    TODO: - Test it
          - Use also a polynomial

    Parameters
    ----------
    x_size
    y_size
    temperature
    counts

    Returns
    -------

    """

    img = np.zeros(shape=(x_size, y_size)) + counts

    sp = Spextrum().black_body_spectrum(temperature=temperature)

    header = fits.Header({"NAXIS": 2,
                          "NAXIS1": x_size,
                          "NAXIS2": y_size,
             #             "CRPIX1": w // 2,
             #             "CRPIX2": h // 2,
             #             "CRVAL1": 0,
             #             "CRVAL2": 0,
             #             "CDELT1": -1 * plate_scale.to(u.deg).value,
             #             "CDELT2": plate_scale.to(u.deg).value,
             #             "CUNIT1": "DEG",
             #             "CUNIT2": "DEG",
             #             "CTYPE1": 'RA---TAN',
             #             "CTYPE2": 'DEC--TAN',
                          "SPEC_REF": 0})

    hdu = fits.ImageHDU(data=img, header=header)
    src = Source()
    src.spectra = [sp]
    src.fields = [hdu]

    return src





