import numpy as np

from astropy.io import fits
from astropy.utils.data import download_file
from astropy import units as u
from astropy.wcs import WCS

from spextra import Spextrum
from scopesim.source.source_templates import Source

from .. import rc
from ..utils import general_utils as gu

from ..utils.exgal_models import GalaxyBase


def source_from_image(image, sed, pixel_scale, amplitude, filter_curve):
    """
    creates a source from an image (numpy 2D array)

    Parameters
    ----------
    image : np.ndarray
      a 2D numpy array
    sed : basestring
       any sed available in the spextra library or a Spextrum object
    pixel_scale : float, u.Quantity
        pixel scale on the sky (u.arcsec, u.arcmin, u.deg).
    amplitude : float, u.Quantity
       magnitude of the spectra in a filter_curve
    filter_curve : basestring
        a filter in the speXtra database or the SVO

    Returns
    -------

    """
    if isinstance(pixel_scale, u.Quantity) is False:
        pixel_scale = pixel_scale * u.arcsec

    if isinstance(sed, Spextrum) is False:
        sp = Spextrum(sed).scale_to_magnitude(amplitude=amplitude, filter_curve=filter_curve)
    else:
        sp = sed

    w, h = image.shape
    x_0 = w // 2
    y_0 = h // 2

    wcs_dict = dict(NAXIS=2,
                    NAXIS1=w+1,
                    NAXIS2=h+1,
                    CRPIX1=x_0,
                    CRPIX2=y_0,
                    CRVAL1=0,
                    CRVAL2=0,
                    CDELT1=-1 * pixel_scale.to(u.deg).value,
                    CDELT2=pixel_scale.to(u.deg).value,
                    CUNIT1="DEG",
                    CUNIT2="DEG",
                    CTYPE1='RA---TAN',
                    CTYPE2='DEC--TAN')

    wcs = WCS(wcs_dict)

    header = fits.Header(wcs.to_header())
    header.update({"SPEC_REF": 0})

    data = image / np.sum(image)
    hdu = fits.ImageHDU(data=data, header=header)
    src = Source(image_hdu=hdu, spectra=sp)

    return src


def source_from_file(filename, pixel_scale, sed, amplitude, filter_curve, cut=0, ext=1, **kwargs):
    """
    Creates a source from a fits image

    Parameters
    ----------
    filename
    pixel_scale
    sed
    amplitude
    filter_curve
    ext : int
       fits extension
    cut
    **kwargs : arguments passed to fits.open

    Returns
    -------

    """
    image = fits.getdata(filename, ext=ext, **kwargs)
    image[image < cut] = 0
    src = source_from_image(image=image, sed=sed, pixel_scale=pixel_scale, amplitude=amplitude, filter_curve=filter_curve)

    return src


def source_from_cube():
    """
    just placeholder for the necessary function. Possible ways to implement it

    - a list of images of one pixel each associated with a different spextrum. Likely expensive to simulate
    - a list of monocromatic images (slices) with scaled counts according to the each pixel spectrum. Each image is associated
    to a very short spectrum.  Likely expensive to simulate but probably gives less problems finding out the boundaries
    - using some dimensionality reduction find the spaxels with similar SED leading to fewer images to simulate.
    Probably difficult to generalize and implement.

    Returns
    -------

    """
    pass


def source_from_moments():
    """
    Place holder for that function.



    Returns
    -------

    """
    pass