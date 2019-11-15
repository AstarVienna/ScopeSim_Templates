from astropy.table import Table
from synphot import SourceSpectrum

from scopesim_templates.basic import misc, galaxy
from scopesim_templates.rc import Source


class TestSpiralTwoComponent:
    def test_returns_source_object(self):
        src = galaxy.spiral_two_component()
        assert isinstance(src, Source)
        assert len(src.fields) == 2
        assert src.fields[0].data.shape == (1200, 1200)
        assert isinstance(src.spectra[1], SourceSpectrum)

    def test_rescale_flux(self):
        src = galaxy.spiral_two_component()


