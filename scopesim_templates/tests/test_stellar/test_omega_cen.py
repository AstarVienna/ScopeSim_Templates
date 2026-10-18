# -*- coding: utf-8 -*-
"""Tests for the omega Cen template and its time evolution."""

import pytest
import numpy as np
from numpy import testing as npt

from scopesim_templates.rc import Source
from scopesim_templates.utils.general_utils import ab_spectrum
from scopesim_templates.stellar import omega_cen, evolve_omega_cen
from scopesim_templates.stellar.omega_centauri import (
    _omega_cen_catalogue, _ClusterPotential, _photometric_mass,
    TURNOFF_MASS, TURNOFF_ABSMAG)
from scopesim_templates.stellar._orbits import to_jyear



# A small, concentrated cluster keeps the offline tests fast (~1e5 draws)
SMALL = dict(fov=6.0, mass=5e4, r_core=0.3, r_tide=10.0, imbh_mass=2e4)


def _catalogue(**kwargs):
    return _omega_cen_catalogue(**{**SMALL, **kwargs})


def _offline_source(**kwargs):
    """Wrap a catalogue in a Source with a flat spectrum (no downloads)."""
    tbl = _catalogue(**kwargs)
    tbl["ref"] = np.zeros(len(tbl), dtype=int)
    tbl["weight"] = 10 ** (-0.4 * tbl["mag"].value)
    meta = dict(tbl.meta)
    tbl.meta.clear()
    src = Source(spectra=[ab_spectrum()], table=tbl)
    src.meta.update(meta)
    src.meta.update({"epoch": "J2030", "epoch_jyear": to_jyear("J2030")})
    return src


class TestStructure:
    def test_default_potential_matches_literature(self):
        pot = _ClusterPotential()
        # workshop validation: sigma_1D(1 pc) = 20.8 km/s (observed 17-22),
        # IMBH influence radius 2.33 pc for 1e5 Msun
        assert pot.dispersion(1.0) == pytest.approx(20.8, abs=0.2)
        assert pot.r_cusp == pytest.approx(2.33, abs=0.02)
        assert pot.enclosed_mass(pot.r_tide) == pytest.approx(3.55e6)

    def test_photometric_mass_anchors_turnoff(self):
        from scopesim_templates.stellar._population import dwarf_absmag
        phot = _photometric_mass(np.array([0.3, 0.5, TURNOFF_MASS]),
                                 TURNOFF_MASS, TURNOFF_ABSMAG)
        npt.assert_allclose(phot[:2], [0.3, 0.5])  # unevolved stars untouched
        assert dwarf_absmag(phot[2], "Mv") == pytest.approx(TURNOFF_ABSMAG,
                                                            abs=0.05)


class TestCatalogue:
    """Offline: sampling and dynamics only, no spectra."""

    def test_same_seed_identical(self):
        a, b = _catalogue(seed=3), _catalogue(seed=3)
        for col in a.colnames:
            npt.assert_array_equal(a[col], b[col])

    def test_different_seed_differs(self):
        a, b = _catalogue(seed=3), _catalogue(seed=4)
        assert len(a) != len(b) or not np.array_equal(a["x"], b["x"])

    def test_seed_none_not_reproducible(self):
        a, b = _catalogue(seed=None), _catalogue(seed=None)
        assert len(a) != len(b) or not np.array_equal(a["x"], b["x"])

    def test_imbh_mass_changes_only_velocities(self):
        a = _catalogue(seed=5, imbh_mass=0.0)
        b = _catalogue(seed=5, imbh_mass=1e5)
        for col in ("x", "y", "z_pc", "mass", "mag", "phase"):
            npt.assert_array_equal(a[col], b[col])
        assert not np.array_equal(a["vx"], b["vx"])
        # a heavier black hole means faster stars near the centre
        inner = np.hypot(a["x"], a["y"]) < 1.0
        speed = lambda t: np.hypot(t["pm_x"], t["pm_y"])[inner]
        assert np.median(speed(b)) > np.median(speed(a))

    def test_positions_inside_field(self):
        tbl = _catalogue(seed=1)
        half = SMALL["fov"] / 2
        assert np.abs(tbl["x"]).max() < half
        assert np.abs(tbl["y"]).max() < half
        # x/y [arcsec] are the projection of x_pc/y_pc at the distance
        npt.assert_allclose(tbl["x"], tbl["x_pc"] / 5430 * 206264.806,
                            rtol=1e-8)

    def test_population(self):
        tbl = _catalogue(seed=2)
        ms = tbl["phase"] == "MS"
        assert tbl["mass"][ms].max() <= TURNOFF_MASS
        assert (tbl["phase"] == "RGB").sum() > 0
        # giants are the brightest stars
        assert tbl["absmag"][~ms].min() < tbl["absmag"][ms].min()
        assert set(tbl["spec_types"][ms]) == {""}   # filled by omega_cen

    def test_mag_limit(self):
        tbl = _catalogue(seed=2, mag_limit=22)
        assert tbl["mag"].max() <= 22

    def test_no_stars_raises(self):
        with pytest.raises(ValueError, match="No cluster stars"):
            _catalogue(fov=1e-4, seed=1)


class TestEvolution:
    """Offline time evolution of a Source built from the catalogue."""

    def test_same_stars_plausible_motion(self):
        src = _offline_source(seed=8)
        later = evolve_omega_cen(src, "J2040")
        a, b = src.fields[0].field, later.fields[0].field
        for col in ("ref", "mass", "mag", "weight", "phase"):
            npt.assert_array_equal(a[col], b[col])
        # over 10 yr, displacement ~ proper motion * dt
        dx = (np.asarray(b["x"]) - np.asarray(a["x"])) * 1e3   # mas
        npt.assert_allclose(dx, np.asarray(a["pm_x"]) * 10.0,
                            rtol=0.05, atol=0.5)
        assert later.meta["epoch_jyear"] == pytest.approx(2040.0)
        assert later.meta["evolved_from_jyear"] == pytest.approx(2030.0)

    def test_round_trip(self):
        src = _offline_source(seed=8)
        back = evolve_omega_cen(evolve_omega_cen(src, "J2100"), "J2030")
        # Not bit-exact: each leg uses the cluster mass enclosed at its own
        # starting radius. The mismatch must stay below 10 micro-arcsec.
        npt.assert_allclose(back.fields[0].field["x"],
                            src.fields[0].field["x"], atol=1e-5)

    def test_input_not_mutated(self):
        src = _offline_source(seed=8)
        x0 = np.array(src.fields[0].field["x"])
        meta0 = dict(src.meta)
        evolve_omega_cen(src, "J2050")
        npt.assert_array_equal(src.fields[0].field["x"], x0)
        assert dict(src.meta) == meta0

    def test_rejects_foreign_source(self):
        from astropy.table import Table
        tbl = Table({"x": [0.0], "y": [0.0], "ref": [0], "weight": [1.0]})
        with pytest.raises(ValueError, match="omega_cen"):
            evolve_omega_cen(Source(spectra=[ab_spectrum()], table=tbl),
                             "J2040")


@pytest.mark.webtest
class TestOmegaCenSource:

    def test_returns_source(self):
        src = omega_cen(seed=1, mag_limit=27, **SMALL)
        assert isinstance(src, Source)
        tbl = src.fields[0].field
        assert tbl["ref"].max() < len(src.spectra)
        assert set(tbl["ref"]) == set(range(len(src.spectra)))
        npt.assert_allclose(tbl["weight"], 10 ** (-0.4 * tbl["mag"]))
        assert "" not in set(tbl["spec_types"])
        assert src.meta["object"] == "omega_cen"
        assert src.meta["seed"] == 1
        assert src.meta["imbh_mass"] == SMALL["imbh_mass"]
        assert src.meta["ra"] == pytest.approx(201.697, abs=1e-3)

    def test_same_seed_same_source(self):
        a = omega_cen(seed=2, mag_limit=27, **SMALL)
        b = omega_cen(seed=2, mag_limit=27, **SMALL)
        for col in ("x", "y", "ref", "weight"):
            npt.assert_array_equal(a.fields[0].field[col],
                                   b.fields[0].field[col])

    def test_time_argument_equals_evolve(self):
        src = omega_cen(seed=3, mag_limit=27, **SMALL)
        direct = omega_cen(seed=3, mag_limit=27, time="J2045", **SMALL)
        evolved = evolve_omega_cen(src, "J2045")
        npt.assert_allclose(direct.fields[0].field["x"],
                            evolved.fields[0].field["x"])
        assert len(evolved.spectra) == len(src.spectra)
