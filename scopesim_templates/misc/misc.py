import numpy as np
import synphot

from astropy.io import fits

from astropy import units as u
from astropy.wcs import WCS

from spextra import Spextrum
from scopesim.source.source_templates import Source


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


def poorman_cube_source(filename=None, hdu=None, ext=1, pixel_scale=None, amplitude=None, filter_curve=None):
    """
    An source from cube that might work with the current implementation of field_of_view objects.
    It basically divides the cube in monocromatic slices (images). The associated spectra are just delta functions

    """
    if filename is not None:
        header = fits.getheader(filename, ext=ext)
        data = fits.getdata(filename, ext=ext)
    elif hdu is not None:
        header = hdu[ext].header
        data = hdu[ext].data
    else:
        raise ValueError("Please define a datacube")

    try:
        bunit = u.Unit(header["BUNIT"])
    except KeyError as e:
        raise("No BUNIT defined", e)

    if pixel_scale is not None:  # assuming lineal
        try:
            pixel_scale = pixel_scale.to(u.deg).value
        except AttributeError:
            pixel_scale = pixel_scale / 3600

    wcs = WCS(header)

    celwcs = wcs.celestial
#    celwcs.wcs.cd = np.array([[-1*pixel_scale, 0], [0, pixel_scale]])  # needs more thought

    celhdr = celwcs.to_header()

    specwcs = wcs.spectral

    zpix = np.arange(specwcs.spectral.array_shape[0])
    waves = specwcs.pixel_to_world(zpix)   # keeping wavelengths in native coordinates
    zero_spec = np.zeros(shape=zpix.shape) * bunit

    hdus = []
    specs = []

    src = Source()
    for i in zpix:
        imgdata = data[i]
        flux = np.sum(imgdata)
        imgdata = imgdata / flux

        celhdr.update({"SPEC_REF": i})
        hdus.append(fits.ImageHDU(data=imgdata, header=celhdr))

        zero_spec[i] = flux * bunit
        sp = synphot.SourceSpectrum(synphot.Empirical1D, points=waves, lookup_table=zero_spec)
        specs.append(sp)
        zero_spec[i] = 0

        src = src + Source(image_hdu=fits.ImageHDU(data=imgdata, header=celhdr),
                           spectra=sp)


#    src = Source(image_hdu=hdul, spectra=specs)

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