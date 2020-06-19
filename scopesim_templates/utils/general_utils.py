from astropy import units as u
from synphot import SourceSpectrum, Empirical1D, ConstFlux1D


def hdu_to_synphot(hdu):

    wave = hdu.data["wavelength"]
    wave_unit = u.Unit(hdu.header["TUNIT1"])
    flux = hdu.data["flux"]
    flux_unit = u.Unit(hdu.header["TUNIT2"])

    spec = SourceSpectrum(Empirical1D, points=wave*wave_unit,
                          lookup_table=flux*flux_unit)

    return spec


def vega_spectrum(mag=0):
    vega = SourceSpectrum.from_vega(cache=True)
    return vega * 10**(-0.4 * mag)


def st_spectrum(mag=0):
    return SourceSpectrum(ConstFlux1D, amplitude=mag*u.STmag)


def ab_spectrum(mag=0):
    return SourceSpectrum(ConstFlux1D, amplitude=mag*u.ABmag)


def function_call_str(func, args):
    func_str = ".".join([func.__module__, func.__name__])
    args_str = ", ".join(["{}={}".format(key, args[key]) for key in args])
    func_call = "{}({})".format(func_str, args_str)

    return func_call
