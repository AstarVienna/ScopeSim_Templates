# -*- coding: utf-8 -*-

import warnings

import pytest
import numpy as np
from numpy import testing as npt
from astropy import units as u

from scopesim_templates.rc import Source
from scopesim_templates.stellar import globular_cluster
from scopesim_templates.stellar import cluster_utils as cu


# Small-enough draws to keep webtest runtime bounded but with enough
# stars for statistical assertions.
DENSITY = 0.5
FOV = 10.0
DM = 13.6  # ~ 5 kpc
SEED = 1234


def _build(**overrides):
    kwargs = dict(density=DENSITY, fov=FOV, distance_modulus=DM, seed=SEED)
    kwargs.update(overrides)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return globular_cluster(**kwargs)


@pytest.mark.webtest
class TestGlobularCluster:

    def test_returns_source(self):
        src = _build()
        assert isinstance(src, Source)
        assert len(src.fields[0].field) == round(DENSITY * FOV ** 2)

    def test_star_count_matches_density_times_fov_squared(self):
        src = _build(density=1.0, fov=8.0)
        assert len(src.fields[0].field) == 64

    def test_positions_within_fov(self):
        src = _build()
        field = src.fields[0].field
        x = np.asarray(field["x"])
        y = np.asarray(field["y"])
        assert np.all(np.abs(x) <= FOV / 2.0 + 1e-9)
        assert np.all(np.abs(y) <= FOV / 2.0 + 1e-9)

    def test_unique_spectrum_per_star(self):
        src = _build()
        field = src.fields[0].field
        n = len(field)
        assert len(src.spectra) == n
        npt.assert_array_equal(field["ref"], np.arange(n, dtype=int))

    def test_weights_match_apparent_mag(self):
        src = _build()
        field = src.fields[0].field
        masses = np.asarray(field["mass"])
        Mv = np.asarray(cu.mass2absmag(masses), dtype=float)
        expected = 10.0 ** (-0.4 * (Mv + DM))
        npt.assert_allclose(np.asarray(field["weight"]), expected, rtol=1e-6)

    def test_rv_shift_applied(self):
        """Two arbitrary stars with the same spectral type have shifted spectra."""
        src = _build(density=2.0, fov=10.0)  # 200 stars: high chance of dup SpT
        field = src.fields[0].field
        spts = np.asarray(field["spec_types"])
        # Find a SpT shared by >= 2 stars
        unique, counts = np.unique(spts, return_counts=True)
        shared_spt = unique[counts >= 2][0]
        idxs = np.where(spts == shared_spt)[0][:2]
        sp_a, sp_b = src.spectra[int(idxs[0])], src.spectra[int(idxs[1])]
        wave = np.linspace(5000, 6000, 200)  # Angstroms
        flux_a = sp_a(wave)
        flux_b = sp_b(wave)
        assert not np.allclose(flux_a, flux_b, rtol=1e-6, atol=0.0)

    def test_seed_reproducibility(self):
        src1 = _build(seed=99)
        src2 = _build(seed=99)
        f1 = src1.fields[0].field
        f2 = src2.fields[0].field
        npt.assert_array_equal(np.asarray(f1["mass"]), np.asarray(f2["mass"]))
        npt.assert_array_equal(np.asarray(f1["x"]), np.asarray(f2["x"]))
        npt.assert_array_equal(np.asarray(f1["weight"]), np.asarray(f2["weight"]))

    def test_zero_stars_raises(self):
        with pytest.raises(ValueError, match="zero stars"):
            globular_cluster(density=0.0, fov=10.0, distance_modulus=DM)

    def test_salpeter_slope_sanity(self):
        """Rough check: alpha = 2.3 +- 0.3 on the high-mass end of a large draw."""
        src = _build(density=4.0, fov=15.0, seed=7)  # ~900 stars
        masses = np.asarray(src.fields[0].field["mass"])
        # Salpeter alpha on dN/dM; we fit log10(N) vs log10(M) above 1 Msun
        m_hi = masses[masses > 1.0]
        if len(m_hi) < 30:
            pytest.skip("not enough high-mass stars to fit slope reliably")
        # Histogram in log mass
        bins = np.logspace(0.0, np.log10(m_hi.max() * 1.01), 8)
        hist, edges = np.histogram(m_hi, bins=bins)
        keep = hist > 0
        x = 0.5 * (edges[:-1] + edges[1:])[keep]
        y = hist[keep]
        # For log-spaced bins, bin count N ~ M^(1-alpha), so slope = 1 - alpha.
        slope, _ = np.polyfit(np.log10(x), np.log10(y), 1)
        alpha_fit = 1.0 - slope
        assert 1.8 < alpha_fit < 3.0
