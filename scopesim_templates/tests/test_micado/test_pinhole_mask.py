"""Test for flatlamp."""
import numpy as np
from astropy.table import Table
from astropy import units as u
from synphot import SourceSpectrum, Empirical1D

from scopesim_templates.micado.pinhole_masks import pinhole_mask
from scopesim_templates.rc import Source


class TestPinholeMask:
    def test_pinhole_mask_returns_source_object(self):
        """ Example from pinhole_mask docstring"""
        dr = np.arange(-25, 26, 5)      # [arcsec]
        x, y = np.meshgrid(dr, dr)
        x, y = x.flatten(), y.flatten()
        waves = np.arange(0.7, 2.5, 0.001) * u.um

        ph_mask = pinhole_mask(x, y, waves, sum_factor=9001)

        assert isinstance(ph_mask, Source)
        assert isinstance(ph_mask.spectra[0], SourceSpectrum)
        assert isinstance(ph_mask.fields[0], Table)
