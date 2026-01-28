# -*- coding: utf-8 -*-

import pytest
import numpy as np
from numpy import testing as npt
from astropy import units as u
import synphot as sp

from spextra import Passband

from scopesim_templates.stellar import star, stars, star_field, star_grid
from scopesim_templates.rc import Source


def source_eq(source_lhs: Source, source_rhs: Source):
    """hacky way to ensure two source"""
    eq = len(source_lhs.spectra) == len(source_rhs.spectra)
    for spectrum_lhs, spectrum_rhs in zip(source_lhs.spectra.values(),
                                          source_rhs.spectra.values()):
        eq = eq and all(spectrum_lhs.waveset == spectrum_rhs.waveset)
    return eq


class TestStar:
    @pytest.mark.webtest
    def test_star_stars_equivalent(self):
        """verify that star and stars yield the same object and have compatible interfaces"""
        src_from_star = star("Generic/Johnson.V", 0, "A0V", 0, 0)
        src_from_stars = stars("Generic/Johnson.V", [0], ["A0V"], [0], [0])
        src_from_star_unit = star("Generic/Johnson.V", 0*u.mag, "A0V", 0, 0)
        src_from_stars_unit = stars("Generic/Johnson.V", [0*u.mag], ["A0V"], [0], [0])

        _ = source_eq(src_from_star, src_from_stars)
        assert source_eq(src_from_star, src_from_stars)
        assert source_eq(src_from_star, src_from_star_unit)
        assert source_eq(src_from_star, src_from_stars_unit)


@pytest.mark.webtest
class TestStars:
    @pytest.mark.parametrize(
        ("filter_src", "filter_obs", "unit", "expected"),
        [("Generic/Johnson.V", "Generic/Johnson.J", sp.units.PHOTLAM, 193),
         ("Generic/Johnson.V", "Generic/Johnson.I", u.Jansky, 2367),
         ("Paranal/HAWKI.J", "Paranal/HAWKI.J", sp.units.PHOTLAM, 187)]
    )
    def test_returns_correct_photon_counts_when_scaled_to_vega_zero(
            self, filter_src, filter_obs, unit, expected):
        src = stars(filter_src, [0]*u.mag, ["A0V"], [0], [0])
        filt = Passband(filter_obs)
        obs = sp.Observation(src.spectra[0], filt)
        npt.assert_allclose(obs.effstim(unit).value, expected, rtol=.01)

    def test_returns_correct_photon_count_when_scaled_to_jansky(self):
        src = stars("Generic/Johnson.V", [3631]*u.Jansky, ["A0V"], [0], [0])
        filt = Passband("Paranal/HAWKI.J")
        obs = sp.Observation(src.spectra[0], filt)
        phs = obs.effstim(sp.units.PHOTLAM).value * src.fields[0]["weight"][0]
        npt.assert_allclose(phs, 187, rtol=.01)

    def test_only_includes_unique_spectra(self):
        spts = ["A0V"]*3 + ["K1II"]*2
        src = stars("Paranal/HAWKI.J", ([0]*5)*u.mag, spts, [0]*5, [0]*5)
        assert len(src.spectra) == 2
        assert np.max(src.fields[0]["ref"]) == 1
        assert np.all(src.fields[0]["weight"] == 1)


@pytest.mark.webtest
class TestStarField:
    def test_returns_source_object_with_minimal_input(self):
        src = star_field(n=2, mmin=0, mmax=5, width=1)
        assert isinstance(src, Source)


@pytest.mark.webtest
class TestStarGrid:
    def test_returns_source_object_with_correct_weighting_without_unit(self):
        src = star_grid(n=3, mmin=0, mmax=5)
        assert isinstance(src, Source)
        assert src.fields[0]["weight"][2] == 0.01
        assert len(src.spectra) == 1

    def test_returns_source_object_with_correct_weighting_in_jansky(self):
        src = star_grid(n=3, mmin=1*u.Jansky, mmax=3631*u.Jansky)
        assert src.fields[0]["weight"][2] == 3631
