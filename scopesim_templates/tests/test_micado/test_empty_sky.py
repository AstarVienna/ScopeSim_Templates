"""Test for empty_sky."""
from astropy.table import Table
from synphot import SourceSpectrum
from scopesim_templates.micado import empty_sky
from scopesim_templates.rc import Source


class TestEmptySky:
    """Test for empty_sky."""

    def test_empty_sky_returns_source_object(self):
        """Test for empty_sky."""
        sky = empty_sky()
        assert isinstance(sky, Source)
        assert isinstance(sky.spectra[0], SourceSpectrum)
        assert isinstance(sky.fields[0], Table)
        assert sky.fields[0]["ref"][0] == 0
