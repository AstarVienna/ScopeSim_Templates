import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy import units as u
from synphot import SourceSpectrum

from scopesim_templates.rc import Source
from scopesim_templates.extragalactic import galaxies

PLOTS = False


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