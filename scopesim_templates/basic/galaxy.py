from astropy.io import fits
from astropy import units as u
from astropy.utils.data import download_file

from .utils import hdu_to_synphot
from .. import rc


def spiral_two_component(extent=60*u.arcsec, fluxes=[0, 0]*u.mag):

    filename = "spiral_two_component.fits"
    url = rc.__config__["!file.server_url"]
    use_cache = rc.__config__["!file.use_cache"]

    path = download_file(remote_url=url+filename, cache=use_cache)
    hdulist = fits.open(path)
    img_ext = hdulist[0].header["IMG_EXT"]
    spec_ext = hdulist[0].header["SPEC_EXT"]

    src = rc.Source()
    src.fields = hdulist[img_ext:spec_ext]
    src.spectra = [hdu_to_synphot(hdu) for hdu in hdulist[spec_ext:]]

    if not isinstance(extent, u.Quantity):
        extent = u.Quantity(extent, u.arcsec, copy=False)
    extent_deg = extent.to(u.deg).value
    for ii in range(len(src.fields)):
        pass


    return src


