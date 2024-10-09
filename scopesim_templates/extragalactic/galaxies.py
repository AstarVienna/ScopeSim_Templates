# -*- coding: utf-8 -*-
"""TBA."""
import pathlib

import numpy as np
import pyckles

from astropy import units as u
from astropy.io import fits
from astropy.utils import deprecated_renamed_argument
from astropy.wcs import WCS

from spextra import Spextrum

from ..rc import Source, im_plane_utils as ipu, \
    ter_curve_utils as tcu
from ..extragalactic import galaxy_utils as gal_utils
from ..extragalactic.exgal_models import GalaxyBase
from ..misc.misc import source_from_array
from ..utils import general_utils as gu
from ..utils.general_utils import add_function_call_str, hdu_to_synphot


__all__ = ["galaxy",
           "galaxy3d",
           "spiral_two_component",
           "elliptical"]


def _galaxy_setup(pixel_scale, r_eff, extend, **kwargs):
    """Construct GalaxyBase (wrapper), kwargs are passed to that class."""
    r_eff <<= u.arcsec
    r_eff_scaled = r_eff.value / pixel_scale.value

    image_size = 2 * r_eff_scaled * extend  # TODO: Needs unit check
    x_0 = image_size // 2
    y_0 = image_size // 2

    x, y = np.meshgrid(np.arange(image_size),
                       np.arange(image_size))

    gal = GalaxyBase(x=x, y=y, x_0=x_0, y_0=y_0, r_eff=r_eff_scaled,
                     amplitude=1, **kwargs)
    return gal


@deprecated_renamed_argument('plate_scale', 'pixel_scale', '0.1')
@add_function_call_str
def galaxy(sed,           # The SED of the galaxy
           z=0,             # redshift
           amplitude=15,           # magnitude
           filter_curve="g",        # passband
           pixel_scale=0.1,   # the plate scale "/pix
           r_eff=2.5,         # effective radius
           n=4,             # sersic index
           ellip=0.1,         # ellipticity
           theta=0,         # position angle
           extend=3,     # extend in units of r_eff
           ra=0.0,
           dec=0.0):
    """
    Create a source object of a galaxy described by its Sersic index and other parameters.

    This function is ideal for imaging or simple spectroscopy.

    Parameters
    ----------
    sed : str or Spextrum
    z : float
        redshift of the galaxy
    r_eff : float
        effective radius of the galaxy in arcsec, it accepts astropy.units
    amplitude : float
        magnitude or flux of the galaxy, it accepts astropy.units
    filter_curve : str
        name of the filter where the magnitude refer to
    pixel_scale : float
        the scale in arcsec/pixel of the instrument
    n : float
        Sersic index of the galaxy
    ellip : float
        ellipticity of the galaxy
    theta : float
        position angle of the galaxy
    extend : float
        Size of the image in units of r_eff
    ra : float, str
        RA of the source in degree, default 0.0
    dec : float, str
        DEC of the source in degree, default 0.0

    Returns
    -------
    src : scopesim.Source
    """
    pixel_scale <<= u.arcsec
    gal = _galaxy_setup(pixel_scale, r_eff, extend,
                        n=n, ellip=ellip, theta=theta)

    if isinstance(sed, str):
        spec = Spextrum(sed).redshift(z=z)
    elif isinstance(sed, Spextrum.__bases__):
        spec = sed.redshift(z=z)

    src = source_from_array(arr=gal.intensity, sed=spec,
                            pixel_scale=pixel_scale, amplitude=amplitude,
                            filter_curve=filter_curve, ra=ra, dec=dec)
    return src


def _get_masked_subsources(gal, ngrid, scaled_sp, header):
    masks = gal.get_masks(ngrid=ngrid)
    intensity = gal.intensity / gal.intensity.sum()
    velocity = gal.velocity
    dispersion = gal.dispersion

    for i, mask in enumerate(masks):
        # logger.debug("creating subsource %d", i)
        data = mask * intensity  # 0 ... 1

        med_vel = velocity[mask].mean()
        med_sig = dispersion[mask].mean()

        spec = scaled_sp.redshift(vel=med_vel).smooth(sigma=med_sig)

        header.update({"SPEC_REF": i})
        hdu = fits.ImageHDU(data=data, header=header)

        yield Source(image_hdu=hdu, spectra=spec)


@deprecated_renamed_argument('plate_scale', 'pixel_scale', '0.1')
@add_function_call_str
def galaxy3d(sed,           # The SED of the galaxy,
             z=0,             # redshift
             amplitude=15,           # magnitude
             filter_curve="g",        # passband
             pixel_scale=0.1,   # the plate scale "/pix
             r_eff=10,         # effective radius
             n=4,             # sersic index
             ellip=0.1,         # ellipticity
             theta=0,         # position angle
             vmax=100,
             sigma=100,
             extend=2,        # extend in units of r_eff
             ngrid=10):       # griding parameter
    """
    Create a simplified 3D map of a galaxy with flux, rotation velocity and velocity dispersion.

    The maps are binned according to the `ngrid` parameter, higher `ngrid`
    will create finer binned fields, but it may increase the computation time.

    The `ngrid` parameter does not specify the number of bins. A ngrid=10 will
    create around 40 independent regions whilst a `ngrid` of 100 will create
    around 2300 regions.

    This function is ideal for spectroscopy

    Parameters
    ----------
    sed : str or Spextrum
        SED of the galaxy, it can be a string or a Spextrum object.
    z : float
        Redshift of the galaxy.
    amplitude : float, u.Quantity
        Magnitude or flux of the galaxy. The spectrum will be re-escaled to
        this value.
    filter_curve : str
        Name of the filter where the `amplitude` is measured.
    pixel_scale : float
        The scale of the image in arcsec/pixel.
    r_eff : float
        Effective radius of the galaxy in arcsec. It accepts astropy.units.
    n : float
        Sersic index of the galaxy.
    ellip : float
        Ellipticity of the galaxy.
    theta : float
        Position angle of the galaxy.
    vmax : float
        Maximum rotation velocity of the galaxy.
    sigma : float
        Velocity dispersion of the galaxy.
    extend : float
        Size of the image in units of `r_eff`.
    ngrid : int
        Gridding parameter for creating of the galaxy.

    Returns
    -------
    src : scopesim.Source
    """
    vmax <<= u.km / u.s
    sigma <<= u.km / u.s
    amplitude <<= u.ABmag
    pixel_scale <<= u.arcsec

    gal = _galaxy_setup(pixel_scale, r_eff, extend,
                        n=n, ellip=ellip, theta=theta,
                        vmax=vmax, sigma=sigma)

    if isinstance(sed, str):
        shift_sp = Spextrum(sed).redshift(z=z)
        scaled_sp = shift_sp.scale_to_magnitude(amplitude=amplitude,
                                                filter_curve=filter_curve)
    elif isinstance(sed, Spextrum):
        scaled_sp = sed
    else:
        raise TypeError("Parameter 'sed' must be Spextrum, SourceSpectrum or "
                        f"str, found {type(sed)=}.")

    # TODO: don't think we really need to call intensity here
    w, h = gal.intensity.shape
    assert w == gal.x.shape[0], f"{w}, {gal.x.shape}"
    assert h == gal.y.shape[1], f"{h}, {gal.y.shape}"

    wcs_dict = {"NAXIS": 2,
                "NAXIS1": 2 * gal.x_0 + 1,
                "NAXIS2": 2 * gal.y_0 + 1,
                "CRPIX1": (w + 1) / 2,
                "CRPIX2": (h + 1) / 2,
                "CRVAL1": 0,
                "CRVAL2": 0,
                "CDELT1": -1 * pixel_scale.to(u.deg).value,
                "CDELT2": pixel_scale.to(u.deg).value,
                "CUNIT1": "deg",
                "CUNIT2": "deg",
                "CTYPE1": "RA---TAN",
                "CTYPE2": "DEC--TAN"}

    wcs = WCS(wcs_dict)

    header = fits.Header(wcs.to_header())
    header.update({"SPEC_REF": 0})

    return sum(_get_masked_subsources(gal, ngrid, scaled_sp, header),
               start=Source())


@add_function_call_str
def spiral_two_component(extent=60*u.arcsec, fluxes=(0, 0), offset=(0, 0)):
    """
    Create a spiral galaxy using NGC1232L as the template.

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
    # params["object"] = "two component spiral galaxy"

    if isinstance(extent, u.Quantity):
        if extent.unit.physical_type == "angle":
            extent = extent.to(u.deg).value
        else:
            raise ValueError("Physical type of extent must be 'angle', but is:"
                             f" {extent.unit.physical_type}.")
    else:
        extent /= 3600.

    filename = "spiral_two_component.fits.bz2"
    path = pathlib.Path(__file__).parent.absolute() / filename
    with fits.open(path) as hdulist:
        img_ext = hdulist[0].header["IMG_EXT"]
        spec_ext = hdulist[0].header["SPEC_EXT"]

        srcs = []
        for i in range(img_ext, spec_ext):
            w, h = hdulist[i].data.shape
            hdulist[i].header["CRPIX1"] = (w + 1) / 2,
            hdulist[i].header["CRPIX2"] = (h + 1) / 2,
            hdulist[i].header["CRVAL1"] = 0
            hdulist[i].header["CRVAL2"] = 0
            hdulist[i].header["CDELT1"] = extent / w
            hdulist[i].header["CDELT2"] = extent / w
            hdulist[i].header["CUNIT1"] = "deg"
            hdulist[i].header["CUNIT2"] = "deg"
            hdulist[i].header["CTYPE1"] = "RA---TAN"
            hdulist[i].header["CTYPE2"] = "DEC--TAN"
            # hdulist[i].header["SPEC_REF"] = hdulist[i].header["SPEC_EXT"] - spec_ext

            i_spec_ext = spec_ext - img_ext + i
            srcs.append(Source(image_hdu=hdulist[i],
                        spectra=[hdu_to_synphot(hdulist[i_spec_ext])]))
    src = sum(srcs[1:], start=srcs[0])

    # TODO: scale image plane according to fluxes
    # TODO: shift header values according to offset

    # src.meta.update(params)
    # Ensure the number of _meta_dicts is the same as the number of _fields.
    # TODO: check if this still works with the new decorator...
    # src._meta_dicts += [{}] * (len(src.fields) - len(src._meta_dicts))

    return src


@deprecated_renamed_argument('magnitude', 'amplitude', '0.1')
@deprecated_renamed_argument('half_light_radius', 'r_eff', '0.1')
# @add_function_call_str
def elliptical(r_eff, pixel_scale, filter_name, amplitude,
               spectrum="NGC_0584", **kwargs):
    """
    Create a extended :class:`.Source` object for an elliptical galaxy.

    .. note:: This docstring is from simcado, needs to be updated

    Parameters
    ----------
    r_eff : float
        [arcsec]

    pixel_scale : float
        [arcsec]

    filter_name : str, TransmissionCurve, optional
        Default is "Ks". Values can be either:
        - the name of a SimCADO filter : see optics.get_filter_set()
        - or a TransmissionCurve containing a user-defined filter

    amplitude : float
        [mag, mag/arcsec2]

    spectrum : str, optional
        The spectrum to be associated with the galaxy. Values can either be:
        - the name of a SimCADO SED spectrum : see get_SED_names()
        - an EmissionCurve with a user defined spectrum


    Optional Parameters (passed to ``sersic_profile``)
    --------------------------------------------------
    n : float, optional
        Default = 4. Sersic index
        - n=1 for exponential (spiral),
        - n=4 for de Vaucouleurs (elliptical)

    ellipticity : float
        Default = 0.5

    angle : float
        [deg] Default = 30. Rotation anti-clockwise from the x-axis

    width, height : int
        [arcsec] Dimensions of the image. Default: 512*pixel_scale

    x_offset, y_offset : float
        [arcsec] The distance between the centre of the profile and the centre
        of the image. Default: (dx,dy) = (0,0)

    normalization : str, optional
        ["total", "half-light", "centre"] Where the profile equals unity
        If normalization equals:
        - "total" : [Default] whole image has brightness of ``amplitude`` [mag]
        - "half-light" : the pixels at the half-light radius have a surface
                         brightness of ``magnitude`` [mag/arcsec2]
        - "centre" : the maximum pixels have a surface brightness of
                     ``magnitude`` [mag/arcsec2]


    Returns
    -------
    galaxy_src : Source

    """
    # """
    # 1 make a sersic profile ImageHDU
    # 2 get spectrum from Brown
    # 3 scale spectrum to amplitude, assume ABmag if no unit
    # 4 make Source object
    # """

    # TODO: apply decorator here as well, somehow...
    params = {"n": 4,
              "ellipticity": 0.5,
              "angle": 30,
              "normalization": "total",
              "width": pixel_scale * 512,
              "height": pixel_scale * 512,
              "x_offset": 0,
              "y_offset": 0,
              "redshift": 0,
              "r_eff": r_eff,
              "pixel_scale": pixel_scale,
              "filter_name": filter_name,
              "amplitude": amplitude,
              "spectrum_name": str(spectrum),
              "rescale_spectrum": True}
    params.update(kwargs)
    params["function_call"] = gu.function_call_str(elliptical, kwargs=params)
    params["object"] = "elliptical galaxy"

    # 1 make a sersic profile ImageHDU
    img = gal_utils.sersic_profile(r_eff=r_eff / pixel_scale,    # everything in terms of pixels
                                   n=params["n"],
                                   ellipticity=params["ellipticity"],
                                   angle=params["angle"],
                                   normalization=params["normalization"],
                                   width=params["width"] / pixel_scale,
                                   height=params["height"] / pixel_scale)

    h_w, h_h = 0.5 * params["width"], 0.5 * params["height"]
    xs = (np.array([-h_w, h_w]) + params["x_offset"]) / 3600.
    ys = (np.array([-h_h, h_h]) + params["y_offset"]) / 3600.
    hdr = ipu.header_from_list_of_xy(xs, ys, pixel_scale / 3600.)
    hdu = fits.ImageHDU(data=img, header=hdr)

    # 2 get spectrum from Brown
    if isinstance(spectrum, str):
        brown_lib = pyckles.SpectralLibrary("brown", return_style="synphot")
        spectrum = brown_lib[spectrum]
        spectrum.z = params["redshift"]

    # 3 scale the spectra and get the weights
    if params["rescale_spectrum"]:
        if not isinstance(amplitude, u.Quantity):
            amplitude <<= u.ABmag
        spectrum = tcu.scale_spectrum(spectrum, filter_name, amplitude)

    # 4 make Source object
    src = Source(spectra=[spectrum], image_hdu=hdu)
    src.meta.update(params)

    return src
