import numpy as np

from astropy.io import fits
import astropy.units as u


def _basic_imagehdu(n=11):
    hdr_dict = {"CRPIX1": n / 2,
                "CRPIX2": n / 2,
                "CRVAL1": 0,
                "CRVAL2": 0,
                "CDELT1": 0.1 / 3600,
                "CDELT2": 0.1 / 3600,
                "CUNIT1":"DEG",
                "CUNIT2":"DEG",
                "CTYPE1": 'RA---TAN',
                "CTYPE2": 'DEC--TAN'}

    im = np.zeros((n,n))
    for x,y in zip([0, 0, -1, -1, 5], [0, -1, 0, -1, 5]):
        im[0] = 1

    hdu = fits.ImageHDU(data=im)
    hdu.header.update(hdr_dict)

    return hdu


def _make_dummy_cube(scale, wave_unit, ref_wave, wave_step, wave_type, bunit):

    if isinstance(scale, u.Quantity) is False:
        scale = scale * u.arcsec
    if isinstance(wave_unit, u.core.Unit) is False:
        wave_unit = u.AA
    if isinstance(ref_wave, u.Quantity) is False:
        ref_wave = ref_wave * u.AA
    if isinstance(wave_step, u.Quantity) is False:
        wave_step = wave_step * u.AA

    data = np.ones(shape=(100, 20, 20))
    header = fits.Header(dict(NAXIS=3,
                              WCSAXES=3,
                              NAXIS1=data.shape[2] + 1,
                              NAXIS2=data.shape[1] + 1,
                              NAXIS3=data.shape[0] + 1,
                              CRPIX1=data.shape[2] // 2,
                              CRPIX2=data.shape[1] // 2,
                              CRPIX3=1,
                              CRVAL1=0,
                              CRVAL2=0,
                              CRVAL3=ref_wave.to(wave_unit).value,
                              CDELT1=-1 * scale.to(u.deg).value,
                              CDELT2=scale.to(u.deg).value,
                              CDELT3=wave_step.to(wave_unit).value,
                              CUNIT1="DEG",
                              CUNIT2="DEG",
                              CUNIT3=wave_unit.to_string(),
                              CTYPE1='RA---TAN',
                              CTYPE2='DEC--TAN',
                              CTYPE3=wave_type,
                              BUNIT=bunit))

    hdu = fits.ImageHDU(data=data, header=header)

    return hdu