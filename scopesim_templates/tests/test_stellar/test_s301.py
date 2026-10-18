# -*- coding: utf-8 -*-
"""Tests for the S301 template (Mang et al. 2026 orbit solutions)."""

import pytest
import numpy as np
from numpy import testing as npt
from astropy import units as u
from astropy.time import Time

from scopesim_templates.rc import Source
from scopesim_templates.stellar import s301, s301_orbit
from scopesim_templates.stellar.s301_star import (BRACKETT_GAMMA,
                                                  S301_SOLUTIONS)

from ._rng_guard import global_rng_untouched  # noqa: F401


C_KMS = 299792.458
PERI_2031 = 2031.806


class TestS301Orbit:
    """Pure-math checks against the numbers quoted in the paper."""

    def test_two_solutions(self):
        assert s301_orbit(1) is S301_SOLUTIONS[1]
        assert s301_orbit(2) is S301_SOLUTIONS[2]

    @pytest.mark.parametrize("bad", [0, 3, "x", None])
    def test_invalid_solution_raises(self, bad):
        with pytest.raises(ValueError):
            s301_orbit(bad)

    def test_semi_major_axis_and_kepler_period(self):
        orbit = s301_orbit(1)
        assert orbit.a_au == pytest.approx(687, abs=1)
        # fitted 8.68 +/- 0.11 yr; Kepler III from a and M_BH gives 8.69
        assert orbit.kepler_period == pytest.approx(8.69, abs=0.01)

    def test_pericentre_distance_and_speed(self):
        # paper: 136-142 R_s, ~25600 / ~25000 km/s
        assert s301_orbit(1).r_peri_rs == pytest.approx(136, abs=1)
        assert s301_orbit(2).r_peri_rs == pytest.approx(142.5, abs=1)
        assert s301_orbit(1).v_peri == pytest.approx(25600, rel=0.01)
        assert s301_orbit(2).v_peri == pytest.approx(25000, rel=0.01)

    def test_next_pericentre_is_late_2031(self):
        for sol in (1, 2):
            peri = s301_orbit(sol).next_pericentre("2026-01-01")
            assert peri == pytest.approx(PERI_2031, abs=0.002)

    def test_rv_sign_flips_between_solutions(self):
        """The two solutions are near-mirrors: opposite Newtonian RVs."""
        # (avoid the ~day around pericentre where the RV swings through 0
        # at slightly different times for the two fits)
        epochs = PERI_2031 + np.array([-300, -30, 1, 2, 30, 300]) / 365.25
        rv1 = s301_orbit(1).radial_velocity(epochs)
        rv2 = s301_orbit(2).radial_velocity(epochs)
        assert np.all(np.sign(rv1) == -np.sign(rv2))
        # and the sky tracks agree to within the fitting noise (< 2 mas)
        sep = np.hypot(*(np.array(s301_orbit(1).sky_position(epochs))
                         - np.array(s301_orbit(2).sky_position(epochs))))
        assert sep.max() < 2.0

    def test_stationary_rv_epoch_about_one_day_after_pericentre(self):
        orbit = s301_orbit(1)
        peri = orbit.next_pericentre("2026-01-01")
        t = peri + np.linspace(-1, 4, 20001) / 365.25
        rv = orbit.apparent_radial_velocity(t)
        t_stat = t[np.argmax(np.abs(rv))]
        assert (t_stat - peri) * 365.25 == pytest.approx(1.0, abs=0.3)
        # the line is near-stationary there: < 100 km/s change over 5 hours
        rv_5h = orbit.apparent_radial_velocity(
            t_stat + np.array([-2.5, 2.5]) / 24 / 365.25)
        assert abs(rv_5h[1] - rv_5h[0]) < 100

    def test_relativistic_terms_at_pericentre(self):
        orbit = s301_orbit(1)
        extra = (orbit.apparent_radial_velocity(PERI_2031)
                 - orbit.radial_velocity(PERI_2031))
        assert float(extra[()]) == pytest.approx(2270, rel=0.05)


@pytest.mark.webtest
class TestS301Source:
    EPOCH = Time(PERI_2031 + 1 / 365.25, format="jyear")

    def test_returns_source_with_sgr_a_star(self):
        src = s301(self.EPOCH)
        assert isinstance(src, Source)
        tbl = src.fields[0].field
        assert list(tbl["name"]) == ["S301", "Sgr A*"]
        assert len(src.spectra) == 2
        assert src.meta["object"] == "s301"
        assert src.meta["ra"] == pytest.approx(266.41684)

    def test_without_sgr_a_star(self):
        src = s301(self.EPOCH, with_sgr_a_star=False)
        assert len(src.fields[0].field) == 1

    def test_positions_are_arcsec_offsets_from_sgr_a_star(self):
        orbit = s301_orbit(1)
        east, north = orbit.sky_position(self.EPOCH.jyear)
        tbl = s301(self.EPOCH).fields[0].field
        npt.assert_allclose([tbl["x"][0], tbl["y"][0]],
                            [east * 1e-3, north * 1e-3], rtol=1e-10)
        npt.assert_allclose([tbl["x"][1], tbl["y"][1]], 0.0)
        assert np.hypot(tbl["x"][0], tbl["y"][0]) < 0.2  # arcsec

    def test_on_axis_puts_star_at_origin(self):
        off = s301(self.EPOCH).fields[0].field
        on = s301(self.EPOCH, on_axis=True).fields[0].field
        npt.assert_allclose([on["x"][0], on["y"][0]], 0.0)
        npt.assert_allclose([on["x"][1], on["y"][1]],
                            [-off["x"][0], -off["y"][0]])

    def test_rv_column_matches_orbit(self):
        for sol in (1, 2):
            tbl = s301(self.EPOCH, solution=sol).fields[0].field
            expected = s301_orbit(sol).apparent_radial_velocity(
                self.EPOCH.jyear)
            assert tbl["rv"][0] == pytest.approx(float(expected[()]))

    def test_brackett_gamma_is_shifted_by_km_per_s(self):
        """Guards against the bare-float (m/s) redshift bug."""
        src = s301(self.EPOCH, solution=1, with_sgr_a_star=False)
        rv = src.fields[0].field["rv"][0]
        expected = BRACKETT_GAMMA.to_value(u.AA) * (1 + rv / C_KMS)
        wave = np.linspace(expected - 100, expected + 100, 4001)
        flux = src.spectra[0](wave * u.AA).value
        assert wave[np.argmin(flux)] == pytest.approx(expected, abs=10)
        # far away from the rest wavelength (rv ~ +16000 km/s)
        assert abs(expected - BRACKETT_GAMMA.to_value(u.AA)) > 500
