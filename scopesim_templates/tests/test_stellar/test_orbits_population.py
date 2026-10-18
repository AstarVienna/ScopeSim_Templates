# -*- coding: utf-8 -*-
"""Offline tests of the orbit and IMF helpers behind the science-case templates."""

import warnings

import pytest
import numpy as np
from numpy import testing as npt
from astropy.time import Time

from scopesim_templates.stellar import _orbits as orb
from scopesim_templates.stellar import _population as pop

from ._rng_guard import global_rng_untouched  # noqa: F401


class TestToJyear:
    def test_time_object(self):
        assert orb.to_jyear(Time(2031.5, format="jyear")) == pytest.approx(2031.5)

    def test_jyear_string(self):
        assert orb.to_jyear("J2031.81") == pytest.approx(2031.81)

    def test_iso_string(self):
        assert orb.to_jyear("2000-01-01T12:00:00") == pytest.approx(2000.0,
                                                                    abs=1e-3)

    def test_float_is_mjd(self):
        assert orb.to_jyear(51544.5) == pytest.approx(2000.0, abs=1e-3)

    def test_decimal_year_like_float_warns(self):
        with pytest.warns(UserWarning, match="MJD"):
            orb.to_jyear(2031.8)


class TestSolveKepler:
    @pytest.mark.parametrize("ecc", [0.0, 0.5, 0.9832, 0.999])
    def test_satisfies_keplers_equation(self, ecc):
        mean = np.linspace(-3 * np.pi, 3 * np.pi, 1001)
        ecc_anom = orb.solve_kepler(mean, ecc)
        wrapped = np.mod(mean + np.pi, 2 * np.pi) - np.pi
        npt.assert_allclose(ecc_anom - ecc * np.sin(ecc_anom), wrapped,
                            atol=1e-10)

    def test_vectorised_eccentricities(self):
        ecc = np.linspace(0, 0.99, 50)
        ecc_anom = orb.solve_kepler(np.full(50, 1.0), ecc)
        npt.assert_allclose(ecc_anom - ecc * np.sin(ecc_anom), 1.0,
                            atol=1e-10)


ORBIT = orb.KeplerOrbit(label="test", a=100.0, e=0.7, incl=50.0, node=30.0,
                        argp=120.0, t_peri=2020.0, period=10.0)


class TestKeplerOrbit:
    def test_returns_after_one_period(self):
        t = np.linspace(2021, 2029, 9)
        npt.assert_allclose(ORBIT.sky_position(t),
                            ORBIT.sky_position(t + ORBIT.period), atol=1e-9)
        npt.assert_allclose(ORBIT.radial_velocity(t),
                            ORBIT.radial_velocity(t + ORBIT.period),
                            rtol=1e-9)

    def test_pericentre_distance(self):
        # at pericentre the 3D radius is a(1-e); projected separation <= it
        sep = ORBIT.separation(ORBIT.t_peri)
        assert sep <= ORBIT.a * (1 - ORBIT.e) + 1e-9
        apo = ORBIT.separation(ORBIT.t_peri + ORBIT.period / 2)
        assert apo <= ORBIT.a * (1 + ORBIT.e) + 1e-9
        assert apo > sep

    def test_mirror_same_sky_opposite_rv(self):
        mirror = ORBIT.mirror()
        t = np.linspace(2020, 2030, 21)
        npt.assert_allclose(ORBIT.sky_position(t), mirror.sky_position(t),
                            atol=1e-9)
        npt.assert_allclose(ORBIT.radial_velocity(t),
                            -mirror.radial_velocity(t), atol=1e-6)

    def test_face_on_orbit_has_no_rv(self):
        face_on = orb.KeplerOrbit("f", 100.0, 0.3, 0.0, 10.0, 20.0,
                                  2020.0, 10.0)
        npt.assert_allclose(face_on.radial_velocity(np.linspace(2020, 2030)),
                            0.0, atol=1e-9)

    def test_relativistic_terms_are_redshifts(self):
        t = np.linspace(2020, 2030, 11)
        extra = (ORBIT.apparent_radial_velocity(t)
                 - ORBIT.radial_velocity(t))
        assert np.all(extra > 0)

    def test_next_pericentre(self):
        assert ORBIT.next_pericentre(2035.0) == pytest.approx(2040.0)
        assert ORBIT.next_pericentre("J2035") == pytest.approx(2040.0)


class TestKeplerPropagate:
    MU = 4 * np.pi**2   # au^3/yr^2 for 1 Msun

    def test_circular_orbit_closes(self):
        pos = np.array([[1.0, 0.0, 0.0]])
        vel = np.array([[0.0, 2 * np.pi, 0.0]])  # circular, period 1 yr
        for frac in (0.25, 0.5, 1.0, 7.3):
            new_pos, _ = orb.kepler_propagate(pos, vel, self.MU, frac)
            angle = 2 * np.pi * frac
            npt.assert_allclose(new_pos[0], [np.cos(angle), np.sin(angle), 0],
                                atol=1e-9)

    def test_energy_and_round_trip(self):
        rng = np.random.default_rng(1)
        pos = rng.normal(size=(200, 3))
        vel = rng.normal(scale=6.0, size=(200, 3))  # bound and unbound mixed
        new_pos, new_vel = orb.kepler_propagate(pos, vel, self.MU, 3.7)

        def energy(p, v):
            return 0.5 * np.sum(v**2, 1) - self.MU / np.linalg.norm(p, axis=1)
        npt.assert_allclose(energy(new_pos, new_vel), energy(pos, vel),
                            rtol=1e-8, atol=1e-8)
        back_pos, back_vel = orb.kepler_propagate(new_pos, new_vel,
                                                  self.MU, -3.7)
        npt.assert_allclose(back_pos, pos, atol=1e-7)
        npt.assert_allclose(back_vel, vel, atol=1e-7)

    def test_inputs_not_modified(self):
        pos = np.array([[1.0, 0.0, 0.0]])
        vel = np.array([[0.0, 6.0, 0.0]])
        orb.kepler_propagate(pos, vel, self.MU, 0.3)
        npt.assert_array_equal(pos, [[1.0, 0.0, 0.0]])
        npt.assert_array_equal(vel, [[0.0, 6.0, 0.0]])


class TestBrokenPowerlaw:
    def test_same_seed_same_masses(self):
        a = pop.sample_broken_powerlaw(1000, 0.08, 60,
                                       np.random.default_rng(5))
        b = pop.sample_broken_powerlaw(1000, 0.08, 60,
                                       np.random.default_rng(5))
        c = pop.sample_broken_powerlaw(1000, 0.08, 60,
                                       np.random.default_rng(6))
        npt.assert_array_equal(a, b)
        assert not np.array_equal(a, c)

    def test_within_limits(self):
        m = pop.sample_broken_powerlaw(10000, 0.1, 0.83,
                                       np.random.default_rng(1))
        assert m.min() >= 0.1 and m.max() <= 0.83

    def test_mean_matches_analytic(self):
        m = pop.sample_broken_powerlaw(200_000, 0.01, 0.85,
                                       np.random.default_rng(2))
        assert m.mean() == pytest.approx(
            pop.broken_powerlaw_mean(0.01, 0.85), rel=0.01)

    def test_known_means(self):
        # values quoted in the MAVIS workshop gotchas (scipy quad of xi)
        assert pop.broken_powerlaw_mean(0.01, 60) == pytest.approx(0.3652,
                                                                   abs=2e-4)
        assert pop.broken_powerlaw_mean(0.01, 0.85) == pytest.approx(
            0.1779, abs=2e-4)

    def test_slope_of_high_mass_tail(self):
        m = pop.sample_broken_powerlaw(200_000, 1.0, 50.0,
                                       np.random.default_rng(3))
        # Salpeter: N(>m) ~ m^-1.3
        ratio = (m > 10).sum() / (m > 2).sum()
        expected = (10**-1.3 - 50**-1.3) / (2**-1.3 - 50**-1.3)
        assert ratio == pytest.approx(expected, rel=0.05)

    def test_zero_stars(self):
        assert pop.sample_broken_powerlaw(0, 0.1, 1,
                                          np.random.default_rng()).size == 0

    def test_out_of_range_raises(self):
        with pytest.raises(ValueError):
            pop.broken_powerlaw_mean(0.001, 1.0)


class TestDwarfAbsmag:
    def test_monotonic_on_main_sequence(self):
        mass = np.linspace(0.1, 1.5, 50)
        for band in ("Mv", "M_J", "M_Ks"):
            assert np.all(np.diff(pop.dwarf_absmag(mass, band)) <= 0)

    def test_sun(self):
        assert pop.dwarf_absmag(1.0, "Mv") == pytest.approx(4.86, abs=0.1)
