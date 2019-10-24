import numpy as np

from astropy.io import ascii as ioascii
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.utils.data import download_file
from astropy import units

from synphot import Empirical1D, SpectralElement, SourceSpectrum

from .. import rc


def spiral_two_component(extent=60*u.arcsec, fluxes=[0, 0]*u.mag):

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
    src.spectra = [hdu_to_synphot(hdu) for hdu in hdulist[spec_ext:]]

    return src


def hdu_to_synphot(hdu):

    wave = hdu.data["wavelength"]
    wave_unit = u.Unit(hdu.header["TUNIT1"])
    flux = hdu.data["flux"]
    flux_unit = u.Unit(hdu.header["TUNIT2"])

    spec = SourceSpectrum(Empirical1D, points=wave*wave_unit,
                          lookup_table=flux*flux_unit)

    return spec

