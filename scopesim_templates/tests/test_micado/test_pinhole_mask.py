"""Test for flatlamp."""

import numpy as np
from astropy.table import Table
from astropy import units as u
from matplotlib import pyplot as plt
from synphot import SourceSpectrum

from scopesim_templates.micado.pinhole_masks import pinhole_mask
from scopesim_templates.rc import Source, load_example_optical_train

PLOTS = False


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
        assert isinstance(ph_mask.fields[0].field, Table)

    def test_pinhole_with_basic_instrument(self):
        dr = np.arange(-25, 26, 5)      # [arcsec]
        x, y = np.meshgrid(dr, dr)
        x, y = x.flatten(), y.flatten()
        waves = np.arange(0.7, 2.5, 0.001) * u.um

        src = pinhole_mask(x, y, waves, sum_factor=9001)
        opt = load_example_optical_train()
        opt.observe(src)

        im = opt.image_planes[0].data

        if PLOTS:
            plt.imshow(im)
            plt.show()

        assert np.max(im) > np.average(im)
