from pytest import approx
import numpy as np

import astropy.units as u
import synphot as sp

from scopesim_templates.stellar import star, stars, star_field, star_grid
from scopesim_templates.rc import Source, ter_curve_utils as tcu


def source_eq(source_lhs: Source, source_rhs: Source):
    """hacky way to ensure two source"""
    eq = len(source_lhs.spectra) == len(source_rhs.spectra)
    for spectrum_lhs, spectrum_rhs in zip(source_lhs.spectra.values(),
                                          source_rhs.spectra.values()):
        eq = eq and all(spectrum_lhs.waveset == spectrum_rhs.waveset)
    return eq


class TestStar:
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


class TestStars:
    def test_returns_correct_photon_counts_when_scaled_to_vega_zero(self):
        src = stars("Generic/Johnson.V", [0]*u.mag, ["A0V"], [0], [0])

        filt = tcu.get_filter("Generic/Johnson.J")
        obs = sp.Observation(src.spectra[0], filt)
        assert obs.effstim(sp.units.PHOTLAM).value == approx(193, rel=0.01)

        filt = tcu.get_filter("Generic/Johnson.I")
        obs = sp.Observation(src.spectra[0], filt)
        assert obs.effstim(u.Jansky).value == approx(2367, rel=0.01)

    def test_returns_correct_photon_count_when_initialised_in_hawki_j(self):
        src = stars("Paranal/HAWKI.J", [0]*u.mag, ["A0V"], [0], [0])
        filt = tcu.get_filter("Paranal/HAWKI.J")
        obs = sp.Observation(src.spectra[0], filt)
        assert obs.effstim(sp.units.PHOTLAM).value == approx(210, rel=0.01)

    def test_returns_correct_photon_count_when_scaled_to_jansky(self):
        src = stars("Generic/Johnson.V", [3631]*u.Jansky, ["A0V"], [0], [0])
        filt = tcu.get_filter("Paranal/HAWKI.J")
        obs = sp.Observation(src.spectra[0], filt)
        phs = obs.effstim(sp.units.PHOTLAM).value * src.fields[0]["weight"][0]
        assert phs == approx(210, rel=0.01)

    def test_only_includes_unique_spectra(self):
        spts = ["A0V"]*3 + ["K1II"]*2
        src = stars("Paranal/HAWKI.J", ([0]*5)*u.mag, spts, [0]*5, [0]*5)
        assert len(src.spectra) == 2
        assert np.max(src.fields[0]["ref"]) == 1
        assert np.all(src.fields[0]["weight"] == 1)


class TestStarField:
    def test_returns_source_object_with_minimal_input(self):
        src = star_field(n=2, mmin=0, mmax=5, width=1)
        assert isinstance(src, Source)


class TestStarGrid:
    def test_returns_source_object_with_correct_weighting_without_unit(self):
        src = star_grid(n=3, mmin=0, mmax=5)
        assert isinstance(src, Source)
        assert src.fields[0]["weight"][2] == 0.01
        assert len(src.spectra) == 1

    def test_returns_source_object_with_correct_weighting_in_jansky(self):
        src = star_grid(n=3, mmin=1*u.Jansky, mmax=3631*u.Jansky)
        assert src.fields[0]["weight"][2] == 3631
