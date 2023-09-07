import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pytest

from astropy import units as u
from synphot import SourceSpectrum

from scopesim_templates.rc import Source
from scopesim_templates.extragalactic import galaxies

PLOTS = False


class TestGalaxy:
    def test_redshift(self):
        sp0 = galaxies.galaxy("kc96/s0", z=0, amplitude=10*u.ABmag,
                              filter_curve="g", pixel_scale=0.05, r_eff=2.5,
                              n=4, ellip=0.5, theta=45, extend=3).spectra[0]
        sp1 = galaxies.galaxy("kc96/s0", z=1, amplitude=10*u.ABmag,
                              filter_curve="g", pixel_scale=0.05, r_eff=2.5,
                              n=4, ellip=0.5, theta=45, extend=3).spectra[0]
        assert sp0.waveset.max() < sp1.waveset.max()
        kwargs = {"filter_curve": "g", "system_name": "AB"}
        assert sp0.get_magnitude(**kwargs) == sp1.get_magnitude(**kwargs)

    @pytest.mark.parametrize("amp", [10*u.ABmag, 15*u.ABmag, 20*u.ABmag])
    def test_scaling(self, amp):
        gal = galaxies.galaxy("kc96/s0", z=0, amplitude=amp,
                              filter_curve="g", pixel_scale=0.05, 
                              r_eff=2.5, n=4, ellip=0.5, theta=45, extend=3)
        kwargs = {"filter_curve": "g", "system_name": "AB"}
        assert gal.spectra[0].get_magnitude(**kwargs) == amp


class TestElliptical:
    def test_it_works(self):
        gal = galaxies.elliptical(125, 1, "V", 10 * u.ABmag, n=1,
                                  spectrum="NGC_0584")
        assert isinstance(gal, Source)

        if PLOTS:
            plt.subplot(121)
            wave = np.arange(0.3, 2.5, 0.001) * u.um
            plt.loglog(wave, gal.spectra[0](wave))

            plt.subplot(122)
            plt.imshow(gal.fields[0].data, norm=LogNorm())
            plt.colorbar()
            plt.show()


class TestSpiralTwoComponent:
    def test_returns_source_object(self):
        src = galaxies.spiral_two_component()
        assert isinstance(src, Source)
        assert len(src.fields) == 2
        assert src.fields[0].data.shape == (1200, 1200)
        assert isinstance(src.spectra[1], SourceSpectrum)

    def test_rescale_flux(self):
        src = galaxies.spiral_two_component()
