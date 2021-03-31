import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.utils.data import download_file
from astropy.utils.decorators import deprecated_renamed_argument

from synphot import Empirical1D, SpectralElement, SourceSpectrum
import pyckles

from .. import rc
from ..rc import ter_curve_utils as tcu
from ..rc import im_plane_utils as ipu

from ..utils.general_utils import function_call_str
from ..utils import general_utils as gen_utils
from ..utils import galaxy_utils as gal_utils


def spiral_two_component(extent=60*u.arcsec, fluxes=(0, 0), offset=(0, 0)):
    """
    Creates a spiral galaxy using NGC1232L as the template

    Two components are included
        - the newer population (spiral arms), and
        - the older popultaion (bulge)

    Parameters
    ----------
    extent : float
        [arcsec]
    fluxes : list of floats, Quantity
        [mag | ABmag | jansky | FLAM | FNU | PHOTLAM | PHOTNU]
    offset : tuple of floats
        [arcsec]

    Returns
    -------
    gal : scopesim.Source

    """
    params = {"extent": extent,
              "fluxes": fluxes,
              "offset": offset}
    pass
    params["function_call"] = gen_utils.function_call_str(spiral_two_component, params)
    params["object"] = "two component spiral galaxy"

    if isinstance(extent, u.Quantity):
        if extent.unit.physical_type == "angle":
            extent = extent.to(u.deg).value
        else:
            raise ValueError("Physical type of extent must be 'angle': "
                             "".format(extent.unit.physical_type))
    else:
        extent /= 3600.

    filename = "spiral_two_component.fits"
    url = rc.__config__["!file.server_url"]
    use_cache = rc.__config__["!file.use_cache"]
    print(url+filename)

    path = download_file(remote_url=url+filename, cache=use_cache)
    hdulist = fits.open(path)
    img_ext = hdulist[0].header["IMG_EXT"]
    spec_ext = hdulist[0].header["SPEC_EXT"]

    src = rc.Source()
    src.fields = hdulist[img_ext:spec_ext]
    src.spectra = [gen_utils.hdu_to_synphot(hdu) for hdu in hdulist[spec_ext:]]

    for ii in range(len(src.fields)):
        w, h = src.fields[ii].data.shape
        src.fields[ii].header["CRPIX1"] = w // 2
        src.fields[ii].header["CRPIX2"] = h // 2
        src.fields[ii].header["CRVAL1"] = 0
        src.fields[ii].header["CRVAL2"] = 0
        src.fields[ii].header["CDELT1"] = extent / w
        src.fields[ii].header["CDELT2"] = extent / w
        src.fields[ii].header["CUNIT1"] = "DEG"
        src.fields[ii].header["CUNIT2"] = "DEG"
        src.fields[ii].header["CTYPE1"] = "RA---TAN"
        src.fields[ii].header["CTYPE2"] = "DEC--TAN"
        src.fields[ii].header["SPEC_REF"] = src.fields[ii].header["SPEC_EXT"] - spec_ext

    # ..todo: scale image plane according to fluxes
    # ..todo: shift header values according to offset

    src.meta.update(params)

    return src


@deprecated_renamed_argument('magnitude', 'amplitude', '0.1')
def elliptical(half_light_radius, pixel_scale, filter_name, amplitude,
               spectrum="NGC_0584", **kwargs):
    """
    Create a extended :class:`.Source` object for an elliptical galaxy

    .. note:: This docstring is from simcado, needs to be updated

    Parameters
    ----------
    half_light_radius : float
        [arcsec]

    pixel_scale : float
        [arcsec]

    amplitude : float
        [mag, mag/arcsec2]

    n : float, optional
        Power law index. Default = 4
        - n=1 for exponential (spiral),
        - n=4 for de Vaucouleurs (elliptical)

    filter_name : str, TransmissionCurve, optional
        Default is "Ks". Values can be either:
        - the name of a SimCADO filter : see optics.get_filter_set()
        - or a TransmissionCurve containing a user-defined filter

    normalization : str, optional
        ["half-light", "centre", "total"] Where the profile equals unity
        If normalization equals:
        - "half-light" : the pixels at the half-light radius have a surface
                         brightness of ``magnitude`` [mag/arcsec2]
        - "centre" : the maximum pixels have a surface brightness of
                     ``magnitude`` [mag/arcsec2]
        - "total" : the whole image has a brightness of ``magnitude`` [mag]

    spectrum : str, optional
        The spectrum to be associated with the galaxy. Values can either be:
        - the name of a SimCADO SED spectrum : see get_SED_names()
        - an EmissionCurve with a user defined spectrum


    Optional Parameters (passed to ``sersic_profile``)
    --------------------------------------------------
    ellipticity : float
        Default = 0.5

    angle : float
        [deg] Default = 30. Rotation anti-clockwise from the x-axis

    width, height : int
        [arcsec] Dimensions of the image. Default: 512*pixel_scale

    x_offset, y_offset : float
        [arcsec] The distance between the centre of the profile and the centre
        of the image. Default: (dx,dy) = (0,0)


    Returns
    -------
    galaxy_src : Source


    See Also
    --------

    """

    # """
    # 1 make a sersic profile ImageHDU
    # 2 get spectrum from Brown
    # 3 scale spectrum to amplitude, assume ABmag if no unit
    # 4 make Source object
    # """

    params = {"n": 4,
              "ellipticity": 0.5,
              "angle": 30,
              "normalization": "total",
              "width": pixel_scale * 512,
              "height": pixel_scale * 512,
              "x_offset": 0,
              "y_offset": 0,
              "redshift": 0,
              "half_light_radius": half_light_radius,
              "pixel_scale": pixel_scale,
              "filter_name": filter_name,
              "amplitude": amplitude,
              "spectrum_name": str(spectrum)}
    params.update(kwargs)
    params["function_call"] = function_call_str(elliptical, params)
    params["object"] = "elliptical galaxy"

    # 1 make a sersic profile ImageHDU
    im = gal_utils.sersic_profile(r_eff=half_light_radius / pixel_scale,    # everything in terms of pixels
                                  n=params["n"],
                                  ellipticity=params["ellipticity"],
                                  angle=params["angle"],
                                  normalization=params["normalization"],
                                  width=params["width"] / pixel_scale,
                                  height=params["height"] / pixel_scale)

    hw, hh = 0.5 * params["width"], 0.5 * params["height"]
    xs = (np.array([-hw, hw]) + params["x_offset"]) / 3600.
    ys = (np.array([-hh, hh]) + params["y_offset"]) / 3600.
    hdr = ipu.header_from_list_of_xy(xs, ys, pixel_scale / 3600.)
    hdu = fits.ImageHDU(data=im, header=hdr)

    # 2 get spectrum from Brown
    if isinstance(spectrum, str):
        brown_lib = pyckles.SpectralLibrary("brown", return_style="synphot")
        spectrum = brown_lib[spectrum]
        spectrum.z = params["redshift"]

    # 3 scale the spectra and get the weights
    if not isinstance(amplitude, u.Quantity):
        amplitude = amplitude << u.ABmag
    spectrum = tcu.scale_spectrum(spectrum, filter_name, amplitude)

    # 4 make Source object
    src = rc.Source(spectra=[spectrum], image_hdu=hdu)
    src.meta.update(params)

    return src

