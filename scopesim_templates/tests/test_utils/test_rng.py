"""Seeding behaviour of the stochastic templates.

The global-RNG guard in ``conftest.py`` additionally checks that none of
these touch NumPy's global random state.
"""

import pytest
import numpy as np
from numpy import testing as npt

from scopesim_templates.rc import __config__
from scopesim_templates.utils.general_utils import get_rng
from scopesim_templates.utils import samplers
from scopesim_templates.stellar import imf, cluster_utils as cu


@pytest.fixture
def config_seed():
    """Temporarily set the package-wide default seed."""
    old = __config__["!random.seed"]
    __config__["!random.seed"] = 9001
    yield 9001
    __config__["!random.seed"] = old


class TestGetRng:
    def test_default_config_seed_is_none(self):
        assert __config__["!random.seed"] is None

    def test_int_seed_is_reproducible(self):
        assert get_rng(42).random() == get_rng(42).random()

    def test_none_is_not_reproducible(self):
        assert get_rng().random() != get_rng().random()

    def test_generator_is_passed_through(self):
        rng = np.random.default_rng(1)
        assert get_rng(rng) is rng

    def test_config_seed_is_used_when_seed_is_none(self, config_seed):
        assert get_rng().random() == get_rng(config_seed).random()


class TestSamplers:
    @staticmethod
    def pdf(x):
        return np.exp(-x)

    @pytest.mark.parametrize("func", [samplers.np_rvs,
                                      samplers.metropolis_sampler])
    def test_seeded_samplers_are_reproducible(self, func):
        npt.assert_array_equal(func(self.pdf, 50, 0, 5, seed=3),
                               func(self.pdf, 50, 0, 5, seed=3))

    def test_randomvariate_is_reproducible_and_counts_trials(self):
        ran1, ntrial = samplers.randomvariate(self.pdf, 50, 0, 5, seed=3)
        ran2, _ = samplers.randomvariate(self.pdf, 50, 0, 5, seed=3)
        npt.assert_array_equal(ran1, ran2)
        # Rejection sampling needs at least as many trials as accepted values
        assert ntrial >= len(ran1) == 50

    def test_np_rvs_returns_n_values_in_range(self):
        ran = samplers.np_rvs(self.pdf, 1000, 0, 5, seed=1)
        assert len(ran) == 1000
        assert ran.min() >= 0 and ran.max() <= 5


class TestImf:
    def test_generate_cluster_seed_is_reproducible(self):
        m1 = imf.Kroupa_2001().generate_cluster(500, seed=7)[0]
        m2 = imf.Kroupa_2001().generate_cluster(500, seed=7)[0]
        npt.assert_array_equal(m1, m2)

    def test_generate_cluster_without_seed_differs(self):
        m1 = imf.Kroupa_2001().generate_cluster(500)[0]
        m2 = imf.Kroupa_2001().generate_cluster(500)[0]
        assert not np.array_equal(m1, m2)

    def test_callers_mass_limits_are_not_mutated(self):
        limits = np.array([0.1, 8.0])
        salpeter = imf.IMF_broken_powerlaw(mass_limits=limits,
                                           powers=np.array([-2.3]))
        salpeter.generate_cluster(5.0, seed=1)   # total < upper mass limit
        npt.assert_array_equal(limits, [0.1, 8.0])


class TestGaussianDistribution:
    def test_seed_is_reproducible(self):
        npt.assert_array_equal(cu.gaussian_distribution(10, 1., seed=5),
                               cu.gaussian_distribution(10, 1., seed=5))


@pytest.mark.webtest
class TestTemplatesSeeding:
    def test_star_field_accepts_seed(self):
        from scopesim_templates.stellar import star_field
        kw = dict(n=5, mmin=10, mmax=15, width=10)
        x1 = star_field(**kw, seed=1).fields[0].field["x"]
        x2 = star_field(**kw, seed=1).fields[0].field["x"]
        x3 = star_field(**kw).fields[0].field["x"]
        npt.assert_array_equal(x1, x2)
        assert not np.array_equal(x1, x3)

    def test_star_field_records_filter_name_in_meta(self):
        from scopesim_templates.stellar import star_field
        src = star_field(n=2, mmin=10, mmax=12, width=5, filter_name="Ks")
        assert src.meta["filter_name"] == "Ks"

    def test_cluster_seed_is_reproducible(self):
        from scopesim_templates.stellar import cluster
        tbl1 = cluster(mass=200, seed=4).fields[0].field
        tbl2 = cluster(mass=200, seed=4).fields[0].field
        npt.assert_array_equal(tbl1["x"], tbl2["x"])
        npt.assert_array_equal(tbl1["weight"], tbl2["weight"])
