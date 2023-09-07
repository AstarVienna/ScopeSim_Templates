from pytest import approx
from collections.abc import Iterable
import numpy as np
from astropy.table import Table
import scopesim_templates.stellar.cluster_utils as cu

from matplotlib import pyplot as plt
PLOTS = False


class TestCatExists:
    def test_cat_parameter_exists(self):
        assert isinstance(cu.MAMAJEK, Table)
        assert len(cu.MAMAJEK) > 0


class TestMass2Mv:
    def test_returns_1_11_for_a0v(self):
        assert cu.mass2Mv(2.3) == approx(1.11, abs=0.1)

    def test_returns_minus_6_for_200_msun(self):
        assert cu.mass2Mv(200) == approx(-6, abs=0.2)

    def test_returns_multiple_mvs_for_masses(self):
        mvs = cu.mass2Mv(np.array([200, 1.02]))
        assert isinstance(mvs, Iterable)
        assert mvs[1] == approx(4.79, rel=0.01)


class TestMass2SpT:
    def test_returns_o5v_for_250_msun(self):
        assert cu.mass2spt(250) == "O5V"

    def test_returns_a0v_for_2_3_msun(self):
        assert cu.mass2spt(2.3) == "A0V"

    def test_returns_g2v_for_1_msun(self):
        assert cu.mass2spt(1.02) == "G2V"

    def test_returns_m95v_for_0_01_msun(self):
        assert cu.mass2spt(0.01) == "M9.5V"

    def test_multiple_masses_at_once(self):
        spts = cu.mass2spt([2.3, 1.02])
        assert spts[1] == "G2V"
        assert isinstance(spts, list)


class TestClosestPickles:
    def test_individual_sources(self):
        assert cu.closest_pickles("O1V") == "O5V"
        assert cu.closest_pickles("G6V") == "G5V"
        assert cu.closest_pickles("G7V") == "G8V"
        assert cu.closest_pickles("M9V") == "M6V"

    def test_multiple_sources(self):
        assert cu.closest_pickles(["O1V", "M9V"])[0] == "O5V"
        assert cu.closest_pickles(["O1V"]*10000)[0] == "O5V"


class TestGaussianDistribution:
    def test_returns_what_we_want(self):
        x, y = cu.gaussian_distribution(10000, 2.35)
        if PLOTS:
            plt.hist(x, bins=40)
            plt.show()

        # Apparently, 0.02 fails occasionally, so use 0.03.
        # TODO: Make these tests deterministic.
        assert np.std(x) == approx(1, rel=0.03)
