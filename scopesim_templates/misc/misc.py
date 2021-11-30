import warnings
import numpy as np
import synphot

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.utils.decorators import deprecated_renamed_argument

from synphot import SourceSpectrum, Empirical1D

from spextra import Spextrum
from ..rc import Source, ter_curve_utils as tu, scopesim_utils as su


def source_from_imagehdu(image_hdu, filter_name, pixel_unit_amplitude=None,
                         inst_pkg_path=None):
    """
    Creates a scopesim.Source object directly from an fits.ImageHDU

    Parameters
    ----------
    image_hdu : fits.ImageHDU
    filter_name : str
        Either a standard filter name or a filter from an instrument package
    pixel_unit_amplitude, optional
        A Quantity that corresponds to a pixel value of 1 in the image
        If not given, header keyword BUNIT is used
    inst_pkg_path : str, optional
        Not yet implemented

    Returns
    -------
    src : scopesim.Source

    Examples
    --------
    Using a generic filter curve from the Spanish VO::

        >>> image_hdu = fits.ImageHDU(data=np.ones((11, 11)))
        >>> # add WCS info to the header here
        >>> filter_name = "Generic/Johnson_UBVRIJHKL.N"
        >>> src = misc.source_from_imagehdu(image_hdu=hdu,
                                            filter_name=filter_name,
                                            pixel_unit_amplitude=20*u.Jy)

    Using the METIS H2O-ice filter from the METIS ScopeSim package::

        >>> import scopesim
        >>> filter_name = scopesim.rc.__search_path__[0] + \
                          "/METIS/filters/TC_filter_H2O-ice.dat"
        >>> src = misc.source_from_imagehdu(image_hdu=hdu,
                                            filter_name=filter_name,
                                            pixel_unit_amplitude=20*u.Jy)

    """
    # if isinstance(inst_pkg_path, str):
    #     import scopesim
    #     scopesim.rc.__search_path__.append(inst_pkg_path)

    if filter_name is None and waverange is None:
        raise ValueError("Wavelength information must be given with either a "
                         "filter_name or a waverange")

    units = image_hdu.header.get("BUNIT")
    amp = 1.

    if isinstance(pixel_unit_amplitude, u.Quantity):
        amp = pixel_unit_amplitude.value
        units = pixel_unit_amplitude.unit
    elif isinstance(pixel_unit_amplitude, (list, np.ndarray)):
        amp = pixel_unit_amplitude

    if units is None:
        raise ValueError("Units must be supplied with either the BUNIT keyword"
                         "or as an astropy.Quantity with pixel_unit_amplitude")

    amp_unit = amp * u.Unit(units)

    if filter_name is not None and isinstance(filter_name, str):
        filter_curve = tu.get_filter(filter_name)
        waves = filter_curve.waverange

    spec = SourceSpectrum(Empirical1D, points=waves, lookup_table=[1, 1])
    tu.scale_spectrum(spec, filter_name=filter_name, amplitude=amp_unit)

    if image_hdu.header.get("SPEC_REF") is None:
        image_hdu.header["SPEC_REF"] = 0

    src = Source(image_hdu=image_hdu, spectra=spec)

    return src


def source_from_imagehdu_with_flux(filename=None, hdu=None, ext=1, pixel_scale=None, flux=None, bunit=None):
    """
    Source from an image where pixel values have units of flux density expressed by bunit.
    It is possible to change the flux and pixel scale to simulate e.g. more distant objects.

    NOTE: This creates a source with a flat spectrum. It is responsibility of the user to use it correctly with
    the proper filter and instrument configuration.
    """

    if filename is not None:
        header = fits.getheader(filename, ext=ext)
        data = fits.getdata(filename, ext=ext)
    elif hdu is not None:
        header = hdu[ext].header
        data = hdu[ext].data
    else:
        raise ValueError("Please define a filename or a ImageHDU")

    if pixel_scale is not None:
        pixel_scale = pixel_scale / 3600
        header.update({"CDELT1": pixel_scale})
        header.update({"CDELT2": pixel_scale})

    if bunit is not None:
        unit = u.Unit(bunit)
    else:
        unit = u.unit(header["BUNIT"])

    if flux is not None and isinstance(flux, u.Quantity) is False:
        total_flux = flux * unit
    elif flux is not None and isinstance(flux, u.Quantity) is True:
        total_flux = flux
    else:
        total_flux = np.sum(data) * unit

    data = data / np.sum(data)

    image_hdu = fits.ImageHDU(data=data, header=header)

    src = Source(mage_hdu=image_hdu, flux=total_flux)

    return src


def source_from_array(arr, sed, pixel_scale, amplitude, filter_curve):
    """
    creates a source from an image (numpy 2D array)
    Parameters
    ----------
    arr : np.ndarray
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

    w, h = arr.shape
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

    data = arr / np.sum(arr)
    hdu = fits.ImageHDU(data=data, header=header)
    src = Source(image_hdu=hdu, spectra=sp)

    return src


@deprecated_renamed_argument(old_name="image", new_name="arr", since="0.3.0",
                             arg_in_kwargs=True)
def source_from_image(**kwargs):
    warnings.warn("source_from_image has been replaced by source_from_array."
                  "This is to avoid confusion with source_from_imagehdu",
                  PendingDeprecationWarning)
    return source_from_array(**kwargs)


source_from_image.__doc__ = source_from_array.__doc__


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