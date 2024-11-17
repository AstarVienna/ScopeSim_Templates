"""Test for darkness."""

from astropy.table import Table
from synphot import SourceSpectrum

from scopesim_templates.micado import darkness
from scopesim_templates.rc import Source


class TestDarkness:
    """Test for darkness."""

    def test_darkness_returns_source_object(self):
        """Test for darkness."""
        sky = darkness()
        assert isinstance(sky, Source)
        assert isinstance(sky.spectra[0], SourceSpectrum)
        assert isinstance(sky.fields[0].field, Table)
        assert sky.fields[0]["ref"][0] == 0
