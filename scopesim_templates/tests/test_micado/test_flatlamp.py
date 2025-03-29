"""Test for flatlamp."""

from astropy.io.fits import ImageHDU
from synphot import SourceSpectrum

from scopesim_templates.micado import flatlamp
from scopesim_templates.rc import Source


class TestFlatlamp:
    """Test for flatlamp."""

    def test_flatlamp_returns_source_object(self):
        """Test for flatlamp."""
        lamp = flatlamp()
        assert isinstance(lamp, Source)
        assert isinstance(lamp.spectra[0], SourceSpectrum)
        assert isinstance(lamp.fields[0].field, ImageHDU)
