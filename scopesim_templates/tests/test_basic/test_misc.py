from astropy.table import Table
from synphot import SourceSpectrum

from scopesim_templates.basic import misc
from scopesim_templates.rc import Source


class TestEmptySky:
    def test_empty_sky_returns_source_object(self):
        sky = misc.empty_sky()
        assert isinstance(sky, Source)
        assert isinstance(sky.spectra[0], SourceSpectrum)
        assert isinstance(sky.fields[0], Table)
        assert sky.fields[0]["ref"][0] == 0


class TestFlatField:
    def test_flatfield_returns_source_object(self):
        flat = misc.flat_field(temperature=1111)
        assert isinstance(flat, Source)
        assert isinstance(flat.spectra[0], SourceSpectrum)
        assert isinstance(flat.fields[0], Table)
        assert flat.fields[0]["ref"][0] == 0

