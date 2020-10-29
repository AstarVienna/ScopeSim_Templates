import pytest
from pytest import approx

import numpy as np
from astropy import units as u
import synphot as sp

from scopesim_templates.basic.galaxy import *
from scopesim_templates.rc import ter_curve_utils as tcu
from scopesim_templates.rc import Source

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False


class TestElliptical:
    def test_it_works(self):
        gal = elliptical(125, 1, "V", 10 * u.ABmag, n=1, spectrum="NGC_0584")
        assert isinstance(gal, rc.Source)

        if PLOTS:
            plt.subplot(121)
            wave = np.arange(0.3, 2.5, 0.001) * u.um
            plt.loglog(wave, gal.spectra[0](wave))

            plt.subplot(122)
            plt.imshow(gal.fields[0].data, norm=LogNorm())
            plt.colorbar()
            plt.show()

