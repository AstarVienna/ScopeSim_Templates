"""Test for flatlamp."""
from astropy.table import Table
from synphot import SourceSpectrum
from scopesim_templates.micado import flatlamp
from scopesim_templates.rc import Source


class TestFlatlamp:
    """Test for flatlamp."""

    def test_flatlamp_returns_source_object(self):
        """Test for flatlamp."""
        sky = flatlamp()
        assert isinstance(sky, Source)
        assert isinstance(sky.spectra[0], SourceSpectrum)
        assert isinstance(sky.fields[0], Table)
        assert sky.fields[0]["ref"][0] == 0
