from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.utils.data import download_file
from astropy import units as u

from synphot import Empirical1D, SpectralElement, SourceSpectrum
import pyckles

from .. import rc
from ..utils import general_utils as gu


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

    return src
