import numpy as np

from astropy.io import fits
from astropy.utils.data import download_file
from astropy import units as u

from spextra import Spextrum
from scopesim.source.source_templates import Source

from .. import rc
from ..utils import general_utils as gu

from .exgal_models import GalaxyBase


def galaxy(sed,           # The SED of the galaxy
           z=0,             # redshift
           mag=15,           # magnitude
           filter_name="g",        # passband
           plate_scale=0.1,   # the plate scale "/pix
           r_eff=2.5,         # effective radius
           n=4,             # sersic index
           ellip=0.1,         # ellipticity
           theta=0,         # position angle
           extend=2):        # extend in units of r_eff

    """
    Creates a source object of a galaxy described by its Sersic index and other
    parameters.

    This function is ideal for imaging

    Parameters
    ----------
    sed : str or Spextrum
    z : float
        redshift of the galaxy
    r_eff : float
        effective radius of the galaxy in arcsec, it accepts astropy.units
    mag : float
        magnitude of the galaxy, it accepts astropy.units
    filter_name : str
        name of the filter where the magnitude refer to
    plate_scale : float
        the scale in arcsec/pixel of the instrument
    n : float
        Sersic index of the galaxy
    ellip : float
        ellipticity of the galaxy
    theta : float
        position angle of the galaxy
    extend : float
        Size of the image in units of r_eff


    Returns
    -------
    src : scopesim.Source
    """

    if isinstance(mag, u.Quantity) is False:
        mag = mag * u.ABmag
    if isinstance(plate_scale, u.Quantity) is False:
        plate_scale = plate_scale * u.arcsec
    if isinstance(r_eff, u.Quantity) is False:
        r_eff = r_eff * u.arcsec
    if isinstance(sed, str):
        sp = Spextrum(sed).redshift(z=z)
        scaled_sp = sp.scale_to_magnitude(amplitude=mag, filter_curve=filter_name)
    elif isinstance(sed, (Spextrum)):
        scaled_sp = sed

    r_eff = r_eff.to(u.arcsec)
    plate_scale = plate_scale.to(u.arcsec)

    image_size = 2 * (r_eff.value * extend / plate_scale.value)  # TODO: Needs unit check
    x_0 = image_size // 2
    y_0 = image_size // 2

    x, y = np.meshgrid(np.arange(image_size),
                       np.arange(image_size))

    galaxy = GalaxyBase(x=x, y=y, x_0=x_0, y_0=y_0,
                        r_eff=r_eff.value, amplitude=1,  n=n,
                        ellip=ellip, theta=theta)

    img = galaxy.intensity

    w, h = img.shape
    header = fits.Header({"NAXIS": 2,
                          "NAXIS1": 2*x_0 + 1,
                          "NAXIS2": 2*y_0 + 1,
                          "CRPIX1": w // 2,
                          "CRPIX2": h // 2,
                          "CRVAL1": 0,
                          "CRVAL2": 0,
                          "CDELT1": -1*plate_scale.to(u.deg).value,
                          "CDELT2": plate_scale.to(u.deg).value,
                          "CUNIT1": "DEG",
                          "CUNIT2": "DEG",
                          "CTYPE1": 'RA---TAN',
                          "CTYPE2": 'DEC--TAN',
                          "SPEC_REF": 0})

    hdu = fits.ImageHDU(data=img, header=header)
    # hdu.writeto("deleteme.fits", overwrite=True)
    src = Source()
    src.spectra = [scaled_sp]
    src.fields = [hdu]

    return src


def galaxy3d(sed,           # The SED of the galaxy,
             z=0,             # redshift
             mag=15,           # magnitude
             filter_name="g",        # passband
             plate_scale=0.2,   # the plate scale "/pix
             r_eff=10,         # effective radius
             n=4,             # sersic index
             ellip=0.1,         # ellipticity
             theta=0,         # position angle
             vmax=100,
             sigma=100,
             extend=2,        # extend in units of r_eff
             ngrid=10):       # griding parameter

    """
    Creates a source object of a galaxy described by its Sersic index and other
    parameters. It also generates a velocity field (set by vmax) and
    a velocity dispersion map (set by sigma).

    The maps are binned according to the ngrid parameter, higher ngrid will create
    finer binned fields but it may increase the computation time.

    The ngrid parameter does not specify the number of bins. A ngrid=10 will create
    around 40 independent regions whilst a ngrid of 100 will create around 2300 regions

    This function is ideal for spectroscopy

    Parameters
    ----------
    sed : str or Spextrum
        SED of the galaxy, it can be a string or a Spextrum object, in the later case it won't
        be re-escaled.
    z : float
        redshift of the galaxy
    mag : float
        magnitude of the galaxy. The spectrum will be re-escaled to this magnitude
    filter_name : str
        name of the filter where the magnitude is measured
    plate_scale : float
        the scale of the image in arcsec/pixel
    r_eff : float
        effective radius of the galaxy in arcsec. It accepts astropy.units
    n : float
        Sersic index of the galaxy
    ellip : float
        ellipticity of the galaxy
    theta : float
        position angle of the galaxy
    vmax : float
        maximum rotation velocity of the galaxy
    sigma : float
        velocity dispersion of the galaxy

    extend : float
        Size of the image in units of r_eff

    ngrid : int
        gridding parameter for creating of the galaxy

    Returns
    -------
    src : scopesim.Source

    """

    if isinstance(mag, u.Quantity) is False:
        mag = mag * u.ABmag
    if isinstance(plate_scale, u.Quantity) is False:
        plate_scale = plate_scale * u.arcsec
    if isinstance(r_eff, u.Quantity) is False:
        r_eff = r_eff * u.arcsec
    if isinstance(vmax, u.Quantity) is False:
        vmax = vmax*u.km/u.s
    if isinstance(sigma, u.Quantity) is False:
        sigma = sigma*u.km/u.s
    if isinstance(sed, str):
        sp = Spextrum(sed).redshift(z=z)
        scaled_sp = sp.scale_to_magnitude(amplitude=mag, filter_curve=filter_name)
    elif isinstance(sed, Spextrum):
        scaled_sp = sed

    r_eff = r_eff.to(u.arcsec)
    plate_scale = plate_scale.to(u.arcsec)
    vmax = vmax.to(u.km/u.s)
    sigma = sigma.to(u.km/u.s)

    image_size = 2 * (r_eff.value * extend / plate_scale.value)  # TODO: Needs unit check
    print(image_size, r_eff)
    x_0 = image_size // 2
    y_0 = image_size // 2

    x, y = np.meshgrid(np.arange(image_size),
                       np.arange(image_size))

    galaxy = GalaxyBase(x=x, y=y, x_0=x_0, y_0=y_0,
                        r_eff=r_eff.value/plate_scale.value,
                        amplitude=1, n=n,
                        ellip=ellip, theta=theta, vmax=vmax, sigma=sigma)

    intensity = galaxy.intensity
    velocity = galaxy.velocity.value
    dispersion = galaxy.dispersion.value
    masks = galaxy.get_masks(ngrid=ngrid)
    w, h = intensity.shape
    header = fits.Header({"NAXIS": 2,
                          "NAXIS1": 2 * x_0 + 1,
                          "NAXIS2": 2 * y_0 + 1,
                          "CRPIX1": w // 2,
                          "CRPIX2": h // 2,
                          "CRVAL1": 0,
                          "CRVAL2": 0,
                          "CDELT1": -1*plate_scale.to(u.deg).value,
                          "CDELT2": plate_scale.to(u.deg).value,
                          "CUNIT1": "DEG",
                          "CUNIT2": "DEG",
                          "CTYPE1": 'RA---TAN',
                          "CTYPE2": 'DEC--TAN',
                          "SPEC_REF": 0})

    src = Source()
    src.fields = []
    src.spectra = []
    total_flux = np.sum(intensity)

    hdulist = []

    for i, m in enumerate(masks):
        data = m * intensity
        factor = np.sum(data) / total_flux

        masked_vel = np.ma.array(velocity, mask=m == 0)
        masked_sigma = np.ma.array(dispersion, mask=m == 0)
        med_vel = np.ma.median(masked_vel)
        med_sig = np.ma.median(masked_sigma)

        spec = scaled_sp.redshift(vel=med_vel).smooth(sigma=med_sig) * factor

        header["SPEC_REF"] = i
        hdu = fits.ImageHDU(data=data, header=header)
        hdulist.append(hdu)
        src.spectra.append(spec)

    src.fields = fits.HDUList(hdulist)

    return src





#-------------------------------------------------

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
    params["function_call"] = gu.function_call_str(spiral_two_component, params)
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
    src.spectra = [gu.hdu_to_synphot(hdu) for hdu in hdulist[spec_ext:]]

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
        src.fields[ii].header["SPEC_REF"] = src.fields[ii].header["SPEC_EXT"] \
                                            - spec_ext

    # ..todo: scale image plane according to fluxes
    # ..todo: shift header values according to offset

    src.meta.update(params)

    return src
