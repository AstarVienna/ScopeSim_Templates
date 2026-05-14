# -*- coding: utf-8 -*-

import warnings

import pytest
import numpy as np
from numpy import testing as npt
from astropy import units as u
from astropy.table import QTable
from astropy.time import Time

from scopesim_templates.rc import Source
from scopesim_templates.stellar import galactic_centre
from scopesim_templates.stellar import _gillessen


EPOCH = Time("2024-01-01")


@pytest.mark.webtest
class TestGalacticCentre:

    def test_returns_source(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # hyperbolic-orbit skip
            src = galactic_centre(EPOCH)
        assert isinstance(src, Source)
        field = src.fields[0].field
        assert len(field) >= 30  # 40 catalogue rows minus a few hyperbolic

    def test_spectra_one_per_star(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            src = galactic_centre(EPOCH)
        n = len(src.fields[0].field)
        assert len(src.spectra) == n
        # ref column is the identity mapping into spectra
        npt.assert_array_equal(
            src.fields[0].field["ref"], np.arange(n, dtype=int)
        )

    def test_two_distinct_spectral_types_used(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            src = galactic_centre(EPOCH)
        assert set(src.fields[0].field["spec_types"]) == {"B0V", "M8III"}

    def test_rv_applied_per_star(self):
        """Two stars of the same class but different RV must have shifted spectra."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            src = galactic_centre(EPOCH)

        field = src.fields[0].field
        # Need two early-type stars with appreciably different RVs.
        # We can reconstruct RVs from the gillessen output for the same epoch.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            ref_tbl = _gillessen.stars_at_time(EPOCH)
        early_mask = np.array([s.strip().lower() == "e" for s in ref_tbl["SpT"]])
        early_idx = np.where(early_mask)[0]
        rv = ref_tbl["RV"].to(u.km / u.s).value[early_idx]
        # pick the pair with the largest RV spread
        i, j = early_idx[np.argmax(rv)], early_idx[np.argmin(rv)]
        sp_i, sp_j = src.spectra[int(i)], src.spectra[int(j)]
        # Sample both spectra on a shared wavelength grid; nonzero difference
        # proves the RV shift took effect.
        wave = np.linspace(2.0e4, 2.4e4, 200)  # K-band, Angstroms
        flux_i = sp_i(wave)
        flux_j = sp_j(wave)
        assert not np.allclose(flux_i, flux_j, rtol=1e-6, atol=0.0)

    def test_weights_match_kmag(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            src = galactic_centre(EPOCH)
            ref_tbl = _gillessen.stars_at_time(EPOCH)
        expected = 10 ** (-0.4 * ref_tbl["Kmag"].to(u.mag).value)
        npt.assert_allclose(
            src.fields[0].field["weight"], expected, rtol=1e-6
        )

    def test_rv_column_in_fields(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            src = galactic_centre(EPOCH)
            ref_tbl = _gillessen.stars_at_time(EPOCH)
        field = src.fields[0].field
        assert "rv" in field.colnames
        npt.assert_allclose(
            np.asarray(field["rv"]),
            ref_tbl["RV"].to(u.km / u.s).value,
            rtol=1e-9,
        )

    def test_positions_are_arcsec(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            src = galactic_centre(EPOCH)
        field = src.fields[0].field
        # Positions are in arcsec offsets from Sgr A* — expect |x|,|y| < a few "
        x = np.asarray(field["x"])
        y = np.asarray(field["y"])
        assert np.max(np.abs(x)) < 30.0
        assert np.max(np.abs(y)) < 30.0

    def test_unknown_spt_warns_and_defaults_to_early(self, monkeypatch):
        """Inject a star with an unrecognised SpT and confirm the warning + fallback."""
        # Real run for one row, then mutate.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            real_tbl = _gillessen.stars_at_time(EPOCH)
        fake = QTable(real_tbl[:1])
        fake["SpT"] = ["?"]

        def fake_stars_at_time(*args, **kwargs):
            return fake

        monkeypatch.setattr(_gillessen, "stars_at_time", fake_stars_at_time)

        with pytest.warns(UserWarning, match="unknown SpT"):
            src = galactic_centre(EPOCH)
        assert src.fields[0].field["spec_types"][0] == "B0V"
