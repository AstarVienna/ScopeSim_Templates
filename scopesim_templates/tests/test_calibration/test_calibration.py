import scopesim_templates.calibration
from astropy.table import Table
from scopesim import Source
from synphot import SourceSpectrum


class TestEmptySky:
    def test_empty_sky_returns_source_object(self):
        sky = scopesim_templates.calibration.calibration.empty_sky()
        assert isinstance(sky, Source)
        assert isinstance(sky.spectra[0], SourceSpectrum)
        assert isinstance(sky.fields[0], Table)
        assert sky.fields[0]["ref"][0] == 0