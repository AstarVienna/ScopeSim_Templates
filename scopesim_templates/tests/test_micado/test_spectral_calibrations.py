# This test suite takes long because it creates a 10 million entry spectrum 5x over

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

        wave = np.arange(0.8, 2.5, 1e-6) * u.um
        flux = src.spectra[0](wave)

        if PLOTS:
            plt.plot(wave, flux)
            plt.show()

        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.fields[0].field, fits.ImageHDU)

    def test_returns_flux_scaled_spectrum(self):
        src1 = mic_sc.line_list()
        src2 = mic_sc.line_list(unit_flux=5)

        wave = np.arange(0.81, 2.5, 1e-7) * u.um
        flux1 = src1.spectra[0](wave).value
        flux2 = src2.spectra[0](wave).value

        assert 5 * flux1.sum() == approx(flux2.sum(), rel=1e-4)

    def test_returns_a_smoothed_spectrum(self):
        src1 = mic_sc.line_list(smoothing_fwhm=None)
        src2 = mic_sc.line_list(smoothing_fwhm=0.005)

        wave = np.arange(0.81, 2.5, 1e-7) * u.um
        flux1 = src1.spectra[0](wave).value
        flux2 = src2.spectra[0](wave).value

        assert flux1.max() > flux2.max()
        assert flux1.sum() == approx(flux2.sum(), rel=2e-2)

    # def test_micado_lines(self):
    #
    #     import numpy as np
    #     from astropy import units as u
    #     from matplotlib import pyplot as plt
    #     from matplotlib.colors import LogNorm
    #
    #     import scopesim as sim
    #
    #     sim.rc.__config__["!SIM.file.local_packages_path"] = "F:/Work/irdb"
    #
    #     cmds = sim.UserCommands(use_instrument="MICADO",
    #                             set_modes=["SCAO", "SPEC_3000x20"])
    #     cmds["!OBS.dit"] = 3600  # [s] Exposure time (more or less)
    #     cmds["!OBS.filter_name_fw1"] = "Spec_HK"  # Spec_IJ, Spec_HK
    #     cmds["!OBS.filter_name_fw2"] = "open"
    #     cmds["!SIM.spectral.spectral_bin_width"] = 1e-4
    #
    #     micado = sim.OpticalTrain(cmds)
    #
    #     # Turn off the effects that are not needed (see below for more on this)
    #     micado["skycalc_atmosphere"].include = False
    #     micado["telescope_reflection"].include = False
    #
    #     # Set to False to test if everything is working, then flip to True
    #     USE_FULL_DETECTOR = False
    #     micado["full_detector_array"].include = USE_FULL_DETECTOR
    #     micado["detector_window"].include = not USE_FULL_DETECTOR
    #
    #     src_lines = mic_sc.line_list(unit_flux=3e23, dwave=1e-1, smoothing_fwhm=None)
    #
    #     micado.observe(src_lines)
    #
    #     plt.figure(figsize=(15, 15))
    #     plt.imshow(micado.image_planes[0].data, norm=LogNorm())
    #     plt.colorbar()
    #     plt.show()
    #
    #     fl = micado.fov_manager._fovs_list
    #     for fov in fl[1:3]:
    #         w = fov.spectra[0].model.points[0]
    #         f = fov.spectra[0].model.lookup_table
    #         plt.plot(w, f)
    #     plt.show()
    #
    #     plt.figure(figsize=(15, 15))
    #     plt.imshow(fl[3].cube.data.sum(axis=1), norm=LogNorm())
    #     plt.show()
