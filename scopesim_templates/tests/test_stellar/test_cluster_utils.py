from pytest import approx, mark
from collections.abc import Iterable
import numpy as np
from astropy.table import Table

import scopesim_templates.stellar.cluster_utils as cu

from matplotlib import pyplot as plt

PLOTS = False


class TestCatExists:
    def test_cat_is_table(self):
        assert isinstance(cu.MAMAJEK, Table)

    def test_cat_is_nonempty(self):
        assert len(cu.MAMAJEK) > 0


class TestMass2AbsMag:
    @mark.parametrize("mass, mv, atol",
                      [(2.3, 1.11, 0.1),
                       (200, -6, 0.2)])
    def test_returns_correct_value(self, mass, mv, atol):
        assert cu.mass2absmag(mass) == approx(mv, abs=atol)

    def test_returns_multiple_mvs_for_masses(self):
        mvs = cu.mass2absmag(np.array([200, 1.02]))
        assert isinstance(mvs, Iterable)
        assert mvs[1] == approx(4.79, rel=0.01)


class TestMass2SpT:
    @mark.parametrize("mass, spt",
                      [(250, "O5V"),
                       (2.3, "A0V"),
                       (1.02, "G2V"),
                       (0.01, "M9.5V")])
    def test_returns_correct_value(self, mass, spt):
        assert cu.mass2spt(mass) == spt

    def test_multiple_masses_at_once(self):
        spts = cu.mass2spt([2.3, 1.02])
        assert spts[1] == "G2V"
        assert isinstance(spts, list)


class TestClosestPickles:
    @mark.parametrize("spt_in, spt_out",
                      [("O1V", "O5V"),
                       ("G6V", "G5V"),
                       ("G7V", "G8V"),
                       ("M9V", "M6V"),
                       ("M2.4V", "M25V")])
    def test_individual_sources(self, spt_in, spt_out):
        assert cu.closest_pickles(spt_in) == spt_out

    @mark.parametrize("spts_in, spt_out",
                      [(["O1V", "M9V"], "O5V"),
                       (["O1V"]*10000, "O5V")])
    def test_multiple_sources(self, spts_in, spt_out):
        assert cu.closest_pickles(spts_in)[0] == spt_out


class TestGaussianDistribution:
    def test_returns_what_we_want(self):
        x, y = cu.gaussian_distribution(10000, 2.35)
        if PLOTS:
            plt.hist(x, bins=40)
            plt.show()

        # Apparently, 0.02 fails occasionally, so use 0.03.
        # TODO: Make these tests deterministic.
        assert np.std(x) == approx(1, rel=0.03)
