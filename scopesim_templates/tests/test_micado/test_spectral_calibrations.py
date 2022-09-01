# This test suite takes long because it creates a 10 million entry spectrum 5x over

import pytest
from pytest import approx

import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt

from synphot import SourceSpectrum
from scopesim_templates.micado import spectral_calibrations as mic_sc

PLOTS = False


class TestLineList:
    def test_returns_source_object_with_everything_as_needed(self):
        src = mic_sc.line_list()

        wave = np.arange(0.8, 2.5, 1e-7) * u.um
        flux = src.spectra[0](wave)

        if not PLOTS:
            plt.plot(wave, flux)
            plt.show()

        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.fields[0], fits.ImageHDU)

    def test_returns_flux_scaled_spectrum(self):
        src1 = mic_sc.line_list()
        src2 = mic_sc.line_list(unit_flux=5)

        wave = np.arange(0.81, 2.5, 1e-7) * u.um
        flux1 = src1.spectra[0](wave).value
        flux2 = src2.spectra[0](wave).value

        assert 5 * flux1.sum() == approx(flux2.sum(), rel=1e-4)

    def test_returns_a_smoothed_spectrum(self):
        src1 = mic_sc.line_list(smoothing_fwhm=None)
        src2 = mic_sc.line_list(smoothing_fwhm=5)

        wave = np.arange(0.81, 2.5, 1e-7) * u.um
        flux1 = src1.spectra[0](wave).value
        flux2 = src2.spectra[0](wave).value

        assert flux1.max() > flux2.max()
        assert flux1.sum() == approx(flux2.sum(), rel=1e-4)
