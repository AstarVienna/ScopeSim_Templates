# -*- coding: utf-8 -*-
"""Templates that could be used to simulate more general sources."""

from pathlib import Path
import warnings

import numpy as np
import synphot

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.utils.decorators import deprecated_renamed_argument

from synphot import SourceSpectrum, Empirical1D

from spextra import Spextrum
from ..rc import Source, ter_curve_utils as tu
from ..utils.general_utils import add_function_call_str, make_img_wcs_header,\
    RA0, DEC0


__all__ = ["point_source",
           "uniform_source",
           "source_from_array",
           "source_from_file",
           "source_from_imagehdu",
           "source_from_imagehdu_with_flux",
           "source_from_cube",
           ]


def _make_spectrum(sed):
    if isinstance(sed, (Spextrum, SourceSpectrum)):
        return Spextrum(modelclass=sed)
    try:
        return Spextrum(sed)
    except TypeError as err:
        raise TypeError("Parameter 'sed' must be Spextrum, SourceSpectrum or "
                        f"str, found {type(sed)=}.") from err


def _make_scaled_spectrum(sed, amplitude, filter_curve):
    spec = _make_spectrum(sed)
    if amplitude is None:
        return spec

    # This might be redundant, check spextra method!
    if not isinstance(amplitude, u.Quantity):
        amplitude = amplitude * u.ABmag

    return spec.scale_to_magnitude(amplitude, filter_curve)


@add_function_call_str
def point_source(sed, amplitude=None, filter_curve="V", x=0, y=0,
                 ra=RA0, dec=DEC0) -> Source:
    """
    Create a point source with an arbitrary spectrum.

    This is similar to ``scopesim_templates.stellar.star``.

    Parameters
    ----------
    sed : str or synphot.Source_Spectrum
        Passing a str will try to download a sed from the speXtra database
        alternatively an user manipulated ``synphot.Source_Spectrum`` or
        compatible object can be provided.
    amplitude : float
        flux or magnitude of the object. The SED will be scaled to that value
        if left to `None` (default) no scaling will be performed
    filter_curve : str
        any astronomical filter in the speXtra or the spanish VO database
    x : float
        X-coordinate on the plane in arcsec
    y : float
        Y-coordinate on the plane in arcsec
    ra : float, str
        RA coordinates of the center of the field (not used)
    dec : float, str
        DEC coordinates of the center of the field  (not used)

    Return
    ------
    src : ``scopesim.Source``
    """
    spec = _make_scaled_spectrum(sed, amplitude, filter_curve)
    src = Source(spectra=spec, x=[x], y=[y], ref=[0], weight=[1])
    return src


@add_function_call_str
def uniform_source(sed, amplitude=None, filter_curve="V", extend=60,
                   ra=RA0, dec=DEC0) -> Source:
    """
    Create an extended uniform source with an arbitrary spectrum.

    This function reates an image with extend^2 pixels with size of 1 arcsec^2
    so provided amplitudes.
    are in flux or magnitudes per arcsec^2

    Parameters
    ----------
    sed : synphot or spextra sed
    amplitude : float
        magnitude or flux (PER ARCSEC^2) of the spectrum in the specified
        filter_curve.
    filter_curve : str
        any filter curve
    extend : int
        extension of the field in arcsec, will always produce a square field
    ra : float
        RA of the field center (not used)
    dec : float
        DEC of field center (not used)

    Return
    ------

    src : ``scopesim.Source``
    """
    data = np.ones(shape=(extend, extend))
    header = make_img_wcs_header(ra=ra, dec=dec, pixel_scale=1,
                                 image_size=data.shape)
    hdu = fits.ImageHDU(header=header, data=data)

    spec = _make_scaled_spectrum(sed, amplitude, filter_curve)
    src = Source(spectra=spec, image_hdu=hdu)
    return src


@add_function_call_str
def source_from_imagehdu(image_hdu, filter_name, pixel_unit_amplitude=None,
                         waverange=None, inst_pkg_path=None) -> Source:
    r"""
    Create a ``scopesim.Source`` object directly from a ``fits.ImageHDU``.

    Parameters
    ----------
    image_hdu : fits.ImageHDU
    filter_name : Path | str
        Either a standard filter name or a filter from an instrument package
    waverange: tuple
        wave_min and wave_max of the spectral range
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

        >>> from scopesim_templates.misc import misc
        >>> image_hdu = fits.ImageHDU(
        ...     data=np.ones((11, 11)),
        ...     header=fits.Header({
        ...         "CUNIT1": "deg",
        ...         "CUNIT2": "deg",
        ...         "CDELT1": 0.001,
        ...         "CDELT2": 0.001,
        ...     }))
        >>> # add WCS info to the header here
        >>> filter_name = "Generic/Johnson_UBVRIJHKL.N"
        >>> src = misc.source_from_imagehdu(image_hdu=image_hdu,
        ...                                 filter_name=filter_name,
        ...                                 pixel_unit_amplitude=20*u.Jy)

    Using the METIS H2O-ice filter from the METIS ScopeSim package::

        >>> import scopesim
        >>> _ = scopesim.download_packages(["METIS"])
        >>> filter_name = scopesim.rc.__search_path__[0] / \
        ...               "METIS/filters/TC_filter_H2O-ice.dat"
        >>> src = misc.source_from_imagehdu(image_hdu=image_hdu,
        ...                                 filter_name=filter_name,
        ...                                 pixel_unit_amplitude=20*u.Jy)

    TODO: Check if the image_hdu has WCS
    """
    if filter_name is None and waverange is None:
        raise ValueError("Wavelength information must be given with either a "
                         "filter_name or a waverange")
    if isinstance(filter_name, Path):
        filter_name = str(filter_name)

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
        waverange = filter_curve.waverange

    waves = np.linspace(waverange[0], waverange[1], num=10)
    lookup_table = [1,] * len(waves)
    spec = SourceSpectrum(Empirical1D, points=waves, lookup_table=lookup_table)
    spec = tu.scale_spectrum(spec, filter_name=filter_name, amplitude=amp_unit)

    if image_hdu.header.get("SPEC_REF") is None:
        image_hdu.header["SPEC_REF"] = 0

    src = Source(image_hdu=image_hdu, spectra=spec)
    return src


@add_function_call_str
def source_from_imagehdu_with_flux(image_hdu=None, filename=None, ext=1,
                                   pixel_scale=None, flux=None,
                                   bunit=None) -> Source:
    """
    Create a source from an image where pixel values have known flux and unit.

    The units of the flux density expressed by 'bunit'.

    It is possible to change the flux and pixel scale to simulate e.g. more
    distant objects.

    Parameters
    ----------
    image_hdu : fits.ImageHDU
        ImageHDU instance or
    filename : str
        A fits filename.
    ext : int
        Extension where the data and header is located in the file.
    pixel_scale : float, optional
        The pixel scale in arcsec if not in the headers or a different one is
        needed.
    flux : float, u.Quantity, optional
        The total flux of the image if the user wants to specify a different
        one than in the image.
    bunit : u.Quantity
        The units of flux if BUNIT is not in the headers or a different one is
        needed.

    Returns
    -------
    src : scopesim.Source


    Note
    ----
    Creates a source with a flat spectrum. It is responsibility of the user to
    use it correctly with the proper filter and instrument configuration.

    """
    if image_hdu is not None:
        header = image_hdu.header
        data = image_hdu.data
    elif filename is not None:
        header = fits.getheader(filename, ext=ext)
        data = fits.getdata(filename, ext=ext)
    else:
        raise ValueError("Please provide a filename or an ImageHDU")

    # we need this because ScopeSim doesn't understand sky coordinates
    header.update({"CRVAL1": 0, "CRVAL2": 0})

    if pixel_scale is not None:
        pixel_scale /= 3600
        header.update({"CDELT1": pixel_scale, "CDELT2": pixel_scale})

    if bunit is not None:
        unit = u.Unit(bunit).to_string()
        header.update({"BUNIT": unit})

    if isinstance(flux, u.Quantity):
        unit = flux.unit.to_string()
        header.update({"BUNIT": unit})
        flux = flux.value
    elif flux is None:
        flux = np.sum(data)

    unit = u.Unit(header["BUNIT"])
    total_flux = flux * unit
    data = data / np.sum(data)
    image_hdu = fits.ImageHDU(data=data, header=header)
    src = Source(image_hdu=image_hdu, flux=total_flux)
    return src


@add_function_call_str
def source_from_array(arr, sed, pixel_scale, amplitude, filter_curve,
                      ra=RA0, dec=DEC0) -> Source:
    """
    Create a source from an image (numpy 2D array).

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
    src : ``scopesim.Source``
    """
    pixel_scale <<= u.arcsec

    header = make_img_wcs_header(ra=ra, dec=dec, pixel_scale=pixel_scale,
                                 image_size=arr.shape)
    header.update({"SPEC_REF": 0})

    data = arr / np.sum(arr)
    hdu = fits.ImageHDU(data=data, header=header)

    spec = _make_scaled_spectrum(sed, amplitude, filter_curve)
    src = Source(image_hdu=hdu, spectra=spec)
    return src


@deprecated_renamed_argument(old_name="image", new_name="arr", since="0.3.0",
                             arg_in_kwargs=True)
def source_from_image(**kwargs):
    warnings.warn("source_from_image has been replaced by source_from_array."
                  "This is to avoid confusion with source_from_imagehdu",
                  PendingDeprecationWarning)
    return source_from_array(**kwargs)


source_from_image.__doc__ = source_from_array.__doc__


@add_function_call_str
def source_from_file(filename, pixel_scale, sed, amplitude, filter_curve,
                     cut=0, ext=1, **kwargs) -> Source:
    """
    Create a source from a fits image.

    TODO: Pass the header (or the WCS) to the source object

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
    hdu = fits.open(filename, **kwargs)
    data = hdu[ext].data
#    header = hdu[ext].header
#    wcs = WCS(header)
#    if wcs.has_celestial:
#        ra, dec = wcs.wcs.crval

    data[data < cut] = 0
    src = source_from_array(arr=data, sed=sed, pixel_scale=pixel_scale,
                            amplitude=amplitude, filter_curve=filter_curve)

    return src


def source_from_cube(cube, ext=1) -> Source:
    """
    Wrapper for Source(cube=cube) for completeness.

    In the future more checks will be performed to translate most cubes into a
    suitable format.

    TODO: Set the cube to a different flux level and different pixel scales.
    Redshift the cube

    Parameters
    ----------
    cube : Filename, ImageHDU or HDUList containing the cube
    ext : extension where the data is located

    """
    src = Source(cube=cube, ext=ext)

    return src


def poorman_cube_source(filename=None, hdu=None, ext=1, pixel_scale=None,
                        amplitude=None, filter_curve=None) -> Source:
    """
    Create a very basic source from cube.

    That source might work with the current implementation of ``field_of_view``
    objects.

    It basically divides the cube in monocromatic slices (images).
    The associated spectra are just delta functions.
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
    except KeyError as err:
        raise ValueError("No BUNIT defined") from err

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
    # keeping wavelengths in native coordinates
    waves = specwcs.pixel_to_world(zpix)
    assert len(waves) > 3, "At least 4 wavelength points are required."
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
        spec = synphot.SourceSpectrum(synphot.Empirical1D, points=waves,
                                      lookup_table=zero_spec)
        specs.append(spec)
        zero_spec[i] = 0

        src = src + Source(image_hdu=fits.ImageHDU(data=imgdata,
                                                   header=celhdr),
                           spectra=spec)


#    src = Source(image_hdu=hdul, spectra=specs)

    return src


def source_from_moments():
    """Place holder for that function."""
    raise NotImplementedError()
