# -*- coding: utf-8 -*-
"""Tests for the hypothetical Galactic Centre background population."""

import warnings

import pytest
import numpy as np
from numpy import testing as npt

from scopesim_templates.rc import Source
from scopesim_templates.stellar import (galactic_centre,
                                        galactic_centre_background)
from scopesim_templates.stellar.gc_background import (
    _background_catalogue, GRAVITY_2022_MASS_LIMIT)



EPOCH = "2031-10-22"


class TestBackgroundCatalogue:
    """Offline: sampling and orbits only, no spectra."""

    def test_same_seed_identical(self):
        a = _background_catalogue(EPOCH, seed=7)
        b = _background_catalogue(EPOCH, seed=7)
        for col in a.colnames:
            npt.assert_array_equal(a[col], b[col])

    def test_different_seed_differs(self):
        a = _background_catalogue(EPOCH, seed=7)
        b = _background_catalogue(EPOCH, seed=8)
        assert len(a) != len(b) or not np.array_equal(a["x"], b["x"])

    def test_seed_none_not_reproducible(self):
        a = _background_catalogue(EPOCH, seed=None)
        b = _background_catalogue(EPOCH, seed=None)
        assert len(a) != len(b) or not np.array_equal(a["x"], b["x"])

    def test_generator_seed(self):
        a = _background_catalogue(EPOCH, seed=np.random.default_rng(3))
        b = _background_catalogue(EPOCH, seed=3)
        npt.assert_array_equal(a["x"], b["x"])

    def test_star_count_scales_with_mass(self):
        small = _background_catalogue(EPOCH, total_mass=100, seed=1)
        large = _background_catalogue(EPOCH, total_mass=3000, seed=1)
        assert large.meta["n_sampled"] / small.meta["n_sampled"] == \
            pytest.approx(30, rel=0.01)
        # Kroupa 0.08-60 Msun: mean mass ~0.556 -> 180 stars per 100 Msun
        assert small.meta["n_sampled"] == 180
        assert len(large) / len(small) == pytest.approx(30, rel=0.2)

    def test_gravity2022_alias(self):
        tbl = _background_catalogue(EPOCH, total_mass="gravity2022", seed=1)
        assert tbl.meta["total_mass"] == GRAVITY_2022_MASS_LIMIT == 3000
        with pytest.raises(ValueError):
            _background_catalogue(EPOCH, total_mass="lots", seed=1)

    def test_no_stars_above_turnoff(self):
        tbl = _background_catalogue(EPOCH, total_mass=3000, seed=2)
        assert tbl["mass"].max() <= 1.2
        assert tbl.meta["n_remnant"] > 0

    def test_positions_are_arcsec_offsets(self):
        tbl = _background_catalogue(EPOCH, total_mass=3000, r_max=1.0,
                                    seed=2)
        r = np.hypot(tbl["x"], tbl["y"])
        assert r.max() <= 2.0      # apocentre <= 2 a <= 2 r_max
        assert np.median(r) < 1.0
        assert abs(np.mean(tbl["x"])) < 0.1 and abs(np.mean(tbl["y"])) < 0.1

    def test_faint_population(self):
        tbl = _background_catalogue(EPOCH, total_mass=3000, seed=2)
        # turn-off caps the main sequence: nothing brighter than ~m_K 19
        assert tbl["m_k"].min() > 19.0
        # workshop: median m_K ~22.6 for the stars brighter than m_K = 24
        bright = np.asarray(tbl["m_k"])[np.asarray(tbl["m_k"]) < 24]
        assert 22.0 < np.median(bright) < 23.2

    def test_mag_limit_is_a_subset(self):
        full = _background_catalogue(EPOCH, total_mass=1000, seed=4)
        cut = _background_catalogue(EPOCH, total_mass=1000, mag_limit=23,
                                    seed=4)
        assert cut["m_k"].max() <= 23
        assert set(cut["star_id"]) <= set(full["star_id"])
        npt.assert_array_equal(full["x"][np.isin(full["star_id"],
                                                 cut["star_id"])], cut["x"])

    def test_same_stars_move_with_time(self):
        t0 = _background_catalogue("J2030.0", total_mass=300, seed=5)
        t1 = _background_catalogue("J2030.1", total_mass=300, seed=5)
        npt.assert_array_equal(t0["star_id"], t1["star_id"])
        npt.assert_array_equal(t0["mass"], t1["mass"])
        # displacement over 0.1 yr ~ proper motion * dt, for slow stars
        dx = (np.asarray(t1["x"]) - np.asarray(t0["x"])) * 1e3   # mas
        pm_dx = np.asarray(t0["pm_x"]) * 0.1
        slow = np.hypot(t0["pm_x"], t0["pm_y"]) < 50
        npt.assert_allclose(dx[slow], pm_dx[slow], atol=0.5)
        # typical proper motions of tens of mas/yr (workshop: 12-67)
        assert 5 < np.median(np.hypot(t0["pm_x"], t0["pm_y"])) < 100

    def test_returns_after_one_period(self):
        t0 = _background_catalogue("J2030.0", total_mass=100, seed=6)
        period = np.asarray(t0["period"])
        i = int(np.argmin(period))
        later = _background_catalogue(f"J{2030.0 + period[i]}",
                                      total_mass=100, seed=6)
        assert later["x"][i] == pytest.approx(t0["x"][i], abs=1e-6)
        assert later["y"][i] == pytest.approx(t0["y"][i], abs=1e-6)


@pytest.mark.webtest
class TestBackgroundSource:

    def test_returns_source(self):
        src = galactic_centre_background(EPOCH, total_mass=300, seed=1)
        assert isinstance(src, Source)
        tbl = src.fields[0].field
        assert {"x", "y", "ref", "weight", "rv", "m_k"} <= set(tbl.colnames)
        assert tbl["ref"].max() < len(src.spectra)
        npt.assert_allclose(tbl["weight"], 10 ** (-0.4 * tbl["m_k"]))
        assert src.meta["seed"] == 1
        assert src.meta["total_mass"] == 300
        assert src.meta["object"] == "galactic_centre_background"

    def test_same_seed_same_source(self):
        a = galactic_centre_background(EPOCH, total_mass=200, seed=9)
        b = galactic_centre_background(EPOCH, total_mass=200, seed=9)
        for col in ("x", "y", "ref", "weight"):
            npt.assert_array_equal(a.fields[0].field[col],
                                   b.fields[0].field[col])


@pytest.mark.webtest
class TestGalacticCentreWithBackground:

    def _build(self, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            return galactic_centre(EPOCH, **kwargs)

    def test_default_has_no_background(self):
        src = self._build()
        assert len(src.fields) == 1
        assert "background_mass" not in src.meta

    def test_background_added_as_second_field(self):
        src = self._build(background_mass=500, seed=3,
                          background_kwargs={"mag_limit": 24})
        assert len(src.fields) == 2
        n_s = len(src.fields[0].field)
        bg = src.fields[1].field
        assert src.meta["background_mass"] == 500
        assert src.meta["background_n_stars"] == len(bg)
        assert src.meta["object"] == "galactic_centre"
        assert src.meta["ra"] == pytest.approx(266.41684)
        # background refs point past the S-star spectra
        assert bg["ref"].min() >= n_s
        assert bg["ref"].max() < len(src.spectra)

    def test_background_reproducible(self):
        a = self._build(background_mass=200, seed=11)
        b = self._build(background_mass=200, seed=11)
        npt.assert_array_equal(a.fields[1].field["x"],
                               b.fields[1].field["x"])
