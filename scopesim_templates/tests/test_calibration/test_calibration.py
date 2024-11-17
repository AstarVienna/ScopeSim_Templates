from astropy.io.fits import ImageHDU
from astropy.table import Table
from astropy import units as u
import numpy as np
from synphot import SourceSpectrum

from scopesim_templates import calibration
from scopesim_templates.rc import Source


class TestEmptySky:
    def test_empty_sky_returns_source_object(self):
        sky = calibration.calibration.empty_sky()
        assert isinstance(sky, Source)
        assert isinstance(sky.spectra[0], SourceSpectrum)
        assert isinstance(sky.fields[0].field, Table)
        assert sky.fields[0]["ref"][0] == 0


class TestFlatField:
    def test_flat_field_returns_source_object(self):
        flatfield = calibration.calibration.flat_field()
        assert isinstance(flatfield, Source)
        assert isinstance(flatfield.spectra[0], SourceSpectrum)
        assert isinstance(flatfield.fields[0].field, ImageHDU)


class TestLamp:
    def test_lamp_returns_source_object(self):
        res = 20000  # Resolution of the instrument
        waves = np.arange(1e4, 2.7e4, 100)
        fwhm = 2.6*(waves/res)*np.ones(shape=waves.shape)  # Nyquist sampled
        flux = 1e-10*np.ones(waves.shape)
        lamp = calibration.calibration.lamp(
            waves=waves,
            fwhm=fwhm,
            fluxes=flux,
        )
        assert isinstance(lamp, Source)
        assert isinstance(lamp.spectra[0], SourceSpectrum)
        assert isinstance(lamp.fields[0].field, ImageHDU)
