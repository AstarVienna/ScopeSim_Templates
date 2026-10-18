# -*- coding: utf-8 -*-
"""Tests for extragalactic.deep_field (Hubble Ultra Deep Field template).

Most tests are offline: they use the 3-galaxy toy catalogue, a flat SED and
a box-car bandpass, so no spextra/SVO downloads are needed. The
``webtest``-marked tests use real SEDs/filters and run the field through
ScopeSim's example optical train.
"""

import shutil
import warnings
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table
from synphot import (ConstFlux1D, Empirical1D, Observation, PowerLawFlux1D,
                     SourceSpectrum, SpectralElement)
from synphot.units import PHOTLAM

from scopesim_templates.rc import Source
from scopesim_templates.extragalactic import hubble_ultra_deep_field
from scopesim_templates.extragalactic import deep_field as df
from scopesim_templates.extragalactic._sersic import (
    sersic_enclosed_fraction, sersic_image)


TOY_CATALOGUE = Path(__file__).parent / "data" / "toy_hudf.fits"
TOY_RA, TOY_DEC = 53.1625, -27.7914     # centre of the toy catalogue

# Offline stand-ins for spextra SEDs and SVO filters
FLAT_SED = SourceSpectrum(ConstFlux1D, amplitude=1 * u.ABmag)
RED_SED = SourceSpectrum(PowerLawFlux1D, amplitude=1 * PHOTLAM,
                         x_0=16000 * u.AA, alpha=-2)
BOX_FILTER = SpectralElement(Empirical1D,
                             points=[14400, 14500, 17500, 17600] * u.AA,
                             lookup_table=[0.0, 1.0, 1.0, 0.0])
OFFLINE = dict(sed_early=RED_SED, sed_late=FLAT_SED, clump_sed=FLAT_SED,
               filter_name=BOX_FILTER)


def _toy(**kwargs):
    kwargs.setdefault("catalogue_path", str(TOY_CATALOGUE))
    kwargs.setdefault("ra_center", TOY_RA)
    kwargs.setdefault("dec_center", TOY_DEC)
    for key, value in OFFLINE.items():
        kwargs.setdefault(key, value)
    return hubble_ultra_deep_field(**kwargs)


def _vendored(**kwargs):
    for key, value in OFFLINE.items():
        kwargs.setdefault(key, value)
    return hubble_ultra_deep_field(**kwargs)


def _image_flux(src):
    return np.array([f.data.sum(dtype=float) for f in src.image_fields])


def _total_flux(src):
    total = _image_flux(src).sum()
    for fld in src.table_fields:
        total += np.sum(fld.field["weight"])
    return total


def _write_single(path, pa=45.0, q=0.2, dx=10.0, dy=5.0, mag=18.0,
                  r_eff=2.0):
    """One-galaxy catalogue offset by (dx East, dy North) arcsec."""
    ra = TOY_RA + dx / 3600 / np.cos(np.deg2rad(TOY_DEC))
    dec = TOY_DEC + dy / 3600
    Table({"ra": [ra], "dec": [dec], "z": [0.5], "m_F160W": [mag],
           "r_eff_arcsec": [r_eff], "sersic_n": [1.0], "axis_ratio": [q],
           "position_angle_deg": [pa]}).write(path, overwrite=True)
    return path


def _moments(image):
    """Centroid and second-moment orientation (deg from +x, CCW)."""
    ny, nx = image.shape
    yy, xx = np.mgrid[:ny, :nx]
    w = np.clip(image, 0, None)
    cx, cy = (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    mxx = (w * (xx - cx) ** 2).sum() / w.sum()
    myy = (w * (yy - cy) ** 2).sum() / w.sum()
    mxy = (w * (xx - cx) * (yy - cy)).sum() / w.sum()
    angle = 0.5 * np.degrees(np.arctan2(2 * mxy, mxx - myy))
    return cx, cy, angle


def _axis_diff(a, b):
    """Difference between two axis orientations, modulo 180 deg."""
    return abs((a - b + 90) % 180 - 90)


# ---------------------------------------------------------------------------
# Sersic rendering
# ---------------------------------------------------------------------------

class TestSersicImage:
    @pytest.mark.parametrize("n", [0.5, 1.0, 2.0, 4.0])
    def test_rendered_fraction_matches_analytic(self, n):
        """Well-sampled profile: box sum close to the analytic enclosed
        fraction (slightly above it, the box corners reach sqrt(2)*n_reff)."""
        img, frac = sersic_image(0.5, n, 1.0, 0.0, pixel_scale=0.01,
                                 n_reff=8, normalise="total")
        expected = sersic_enclosed_fraction(n, 8.0)
        assert img.sum() == pytest.approx(frac)
        assert expected - 1e-3 < frac < expected + 0.03
        assert frac <= 1.0 + 1e-3

    def test_sum_normalisation_conserves_flux(self):
        for r_eff, n in [(0.02, 4.0), (0.3, 1.0), (1.0, 6.0)]:
            img, _ = sersic_image(r_eff, n, 0.4, 30.0, pixel_scale=0.06)
            assert img.sum() == pytest.approx(1.0, rel=1e-10)
            assert img.shape[0] % 2 == 1 and img.shape[0] == img.shape[1]

    def test_core_oversampling_converged(self):
        """A cuspy n=4 core: 8x oversampling agrees with 32x."""
        a, _ = sersic_image(0.1, 4.0, 1.0, 0.0, 0.004, oversample=8,
                            normalise="total")
        b, _ = sersic_image(0.1, 4.0, 1.0, 0.0, 0.004, oversample=32,
                            normalise="total")
        assert abs(a.sum() - b.sum()) < 1e-3

    def test_angle_is_ccw_from_x(self):
        img, _ = sersic_image(1.0, 1.0, 0.2, 30.0, pixel_scale=0.05)
        _, _, angle = _moments(img)
        assert _axis_diff(angle, 30.0) < 1.0

    def test_invalid_parameters_raise(self):
        with pytest.raises(ValueError):
            sersic_image(0.0, 1.0, 0.5, 0.0, 0.1)
        with pytest.raises(ValueError):
            sersic_image(1.0, 1.0, 1.5, 0.0, 0.1)


# ---------------------------------------------------------------------------
# Catalogue handling
# ---------------------------------------------------------------------------

class TestCatalogue:
    def test_vendored_catalogue_columns_and_attribution(self):
        cat = df.load_hudf_catalogue()
        assert len(cat) == 3152
        for col in ("x", "y", "m_F160W", "r_eff_arcsec", "sersic_n",
                    "axis_ratio", "position_angle_deg", "log_mstar", "z"):
            assert col in cat.colnames
        assert np.all(cat["z"] > 0)
        assert np.all((cat["axis_ratio"] > 0) & (cat["axis_ratio"] <= 1))
        readme = (Path(df.__file__).parent / "data" / "README.md").read_text(
            encoding="utf-8")
        assert "van der Wel" in readme and "Skelton" in readme
        assert "van der Wel" in df.__doc__ and "Skelton" in df.__doc__

    def test_offsets_x_is_east(self):
        cat = df.load_hudf_catalogue(TOY_CATALOGUE, TOY_RA, TOY_DEC)
        east = cat["ra"] > TOY_RA + 1e-6
        assert np.all(cat["x"][east] > 0)
        np.testing.assert_allclose(cat["x"], [0, 5, -5], atol=1e-3)
        np.testing.assert_allclose(cat["y"], [0, 5, -5], atol=1e-3)

    def test_workshop_column_aliases(self, tmp_path):
        tbl = Table.read(TOY_CATALOGUE)
        tbl.rename_columns(
            ["m_F160W", "r_eff_arcsec", "sersic_n", "axis_ratio",
             "position_angle_deg"], ["mag", "re", "n", "q", "pa"])
        tbl["logm"] = [9.0, 10.0, 11.0]
        path = tmp_path / "workshop.fits"
        tbl.write(path)
        cat = df.load_hudf_catalogue(path)
        assert list(cat["m_F160W"]) == list(Table.read(TOY_CATALOGUE)
                                             ["m_F160W"])
        assert list(cat["log_mstar"]) == [9.0, 10.0, 11.0]

    def test_missing_column_raises(self, tmp_path):
        tbl = Table.read(TOY_CATALOGUE)
        tbl.remove_column("sersic_n")
        path = tmp_path / "bad.fits"
        tbl.write(path)
        with pytest.raises(KeyError, match="sersic_n"):
            df.load_hudf_catalogue(path)

    def test_select_field_does_not_mutate(self):
        cat = df.load_hudf_catalogue(TOY_CATALOGUE, TOY_RA, TOY_DEC)
        x_before = np.array(cat["x"])
        sel = df.select_field(cat, 6.0, x0=5.0, y0=5.0)
        assert len(sel) == 1 and sel["id"][0] == 1
        assert abs(sel["x"][0]) < 1e-3
        np.testing.assert_array_equal(cat["x"], x_before)
        assert len(df.select_field(cat, 6.0, margin=3.0)) == 3
        assert len(df.select_field(cat, None, mag_limit=24)) == 2

    def test_richest_field_matches_brute_force(self):
        cat = df.load_hudf_catalogue()
        size, step = 20.0, 5.0
        x0, y0, score = df.richest_field(cat, size, step=step,
                                         max_offset=40.0)
        bright = cat[cat["m_F160W"] < 25.0]
        grid = np.arange(-40.0, 40.0 + 0.5 * step, step)
        best = -1
        for gx in grid:
            for gy in grid:
                inside = ((np.abs(bright["x"] - gx) < size / 2)
                          & (np.abs(bright["y"] - gy) < size / 2))
                best = max(best, np.sum(10 ** (-0.4 * bright["m_F160W"]
                                               [inside])))
        assert score == pytest.approx(best)
        assert abs(x0) <= 40 and abs(y0) <= 40


# ---------------------------------------------------------------------------
# Clumps and RNG
# ---------------------------------------------------------------------------

class TestClumps:
    @pytest.fixture
    def galaxies(self):
        cat = df.load_hudf_catalogue()
        return df.select_field(cat, 40.0)

    def test_clump_fractions_and_positions(self, galaxies):
        clumps = df.clump_catalogue(galaxies, seed=1)
        taken = clumps.meta["taken"]
        assert len(clumps) > 0
        clumpy = taken > 0
        # rescale=True: each clumpy galaxy gives away 10-25 % exactly
        assert np.all((taken[clumpy] >= 0.10 - 1e-12)
                      & (taken[clumpy] <= 0.25 + 1e-12))
        # only disky galaxies are clumpy
        assert np.all(galaxies["sersic_n"][clumpy] <= 2.5)
        # clumps stay within n_reff * r_eff of their host
        host_re = galaxies["r_eff_arcsec"][clumps["host"]]
        assert np.all(np.hypot(clumps["dx"], clumps["dy"])
                      <= 8.0 * host_re + 1e-9)
        # magnitudes consistent with fractions
        host_mag = galaxies["m_F160W"][clumps["host"]]
        np.testing.assert_allclose(
            clumps["mag"], host_mag - 2.5 * np.log10(clumps["frac"]))

    def test_clumps_follow_position_angle(self):
        """Clumps of a very flat galaxy line up along its major axis."""
        gal = Table({"id": [0], "x": [0.0], "y": [0.0], "z": [2.0],
                     "m_F160W": [20.0], "r_eff_arcsec": [1.0],
                     "sersic_n": [1.0], "axis_ratio": [0.05],
                     "position_angle_deg": [30.0], "log_mstar": [9.0]})
        clumps = df.clump_catalogue(gal, flux_fraction=0.5,
                                    clumpy_table={99: [(0, 1.0)]},
                                    fraction_range=(0.005, 0.01), seed=3)
        angle = np.degrees(np.arctan2(clumps["dy"], clumps["dx"]))
        # PA 30 deg east of north -> 60 deg CCW from +x (East)
        assert np.median(_axis_diff(angle, 60.0)) < 5.0

    def test_seed_determinism(self, galaxies):
        a = df.clump_catalogue(galaxies, seed=42)
        b = df.clump_catalogue(galaxies, seed=42)
        c = df.clump_catalogue(galaxies, seed=43)
        np.testing.assert_array_equal(a["x"], b["x"])
        assert len(a) != len(c) or not np.array_equal(a["x"], c["x"])
        gen = df.clump_catalogue(galaxies,
                                 seed=np.random.default_rng(42))
        np.testing.assert_array_equal(a["frac"], gen["frac"])

    def test_global_rng_untouched(self, galaxies):
        state = np.random.get_state()
        df.clump_catalogue(galaxies, seed=7)
        df.clump_catalogue(galaxies)
        after = np.random.get_state()
        assert state[0] == after[0]
        np.testing.assert_array_equal(state[1], after[1])
        assert state[2:] == after[2:]

    def test_clumpy_fraction_table(self):
        frac = df.clumpy_fraction([9.0, 10.0, 11.0, np.nan],
                                  [2.0, 0.5, 0.5, 2.0])
        np.testing.assert_allclose(frac, [0.60, 0.40, 0.15, 0.60])


# ---------------------------------------------------------------------------
# The Source (offline: flat SEDs, box filter)
# ---------------------------------------------------------------------------

class TestHUDFSource:
    def test_returns_source(self):
        src = _toy()
        assert isinstance(src, Source)
        assert len(src.image_fields) == 3
        assert src.meta["n_galaxies"] == 3
        assert src.meta["object"] == "hubble_ultra_deep_field"
        assert "function_call" in src.meta

    def test_crval_are_offsets_not_absolute(self):
        src = _toy()
        cat = df.load_hudf_catalogue(TOY_CATALOGUE, TOY_RA, TOY_DEC)
        for fld, row in zip(src.image_fields, cat):
            assert fld.header["CRVAL1"] * 3600 == pytest.approx(row["x"],
                                                                abs=1e-6)
            assert fld.header["CRVAL2"] * 3600 == pytest.approx(row["y"],
                                                                abs=1e-6)
            assert abs(fld.header["CRVAL1"]) < 10 / 3600
        assert src.meta["ra_center"] == TOY_RA
        assert src.meta["dec_center"] == TOY_DEC

    def test_vendored_field_offsets_small(self):
        src = _vendored(fov=20.0)
        crval = np.array([[f.header["CRVAL1"], f.header["CRVAL2"]]
                          for f in src.image_fields]) * 3600
        assert np.all(np.abs(crval) <= 10.0 + 1e-6)
        assert src.meta["n_galaxies"] == len(src.image_fields) > 10

    def test_flux_matches_catalogue(self):
        """sum(image) == 10**(-0.4 m) per galaxy (conserve_flux=True)."""
        src = _toy()
        mags = np.array([f.meta["m_F160W"] for f in src.image_fields])
        np.testing.assert_allclose(_image_flux(src), 10 ** (-0.4 * mags),
                                   rtol=1e-5)

    def test_amplitude_passes_through(self):
        src = _toy()
        flux = _image_flux(src)
        mags = Table.read(TOY_CATALOGUE)["m_F160W"]
        assert np.argmax(flux) == np.argmin(mags)
        assert flux.max() / flux.min() == pytest.approx(
            10 ** (0.4 * (mags.max() - mags.min())), rel=1e-5)

    def test_spectra_normalised_to_zero_mag(self):
        src = _toy()
        for spec in src.spectra.values():
            mag = Observation(spec, BOX_FILTER).effstim(u.ABmag).value
            assert mag == pytest.approx(0.0, abs=1e-3)

    def test_conserve_flux_false_keeps_surface_brightness(self):
        src = _toy(conserve_flux=False)
        assert np.all(_image_flux(src) < 10 ** (-0.4 * np.array(
            [f.meta["m_F160W"] for f in src.image_fields])))
        assert 0.9 < src.meta["rendered_fraction_min"] < 1.0

    def test_positions_within_fov(self):
        assert _toy(fov=6.0).meta["n_galaxies"] == 1
        assert _toy(fov=6.0, margin=3.0).meta["n_galaxies"] == 3

    def test_early_late_split_by_n_cut(self):
        """Toy catalogue has n = 4, 1, 3: two early, one late type."""
        src = _toy(n_cut=2.5)
        assert src.meta["n_early"] == 2 and src.meta["n_late"] == 1
        # The two SEDs have different slopes, independent of redshift:
        # flat f_nu -> photlam ~ 1/lam; RED_SED -> photlam ~ lam**2
        expected = {0: (12 / 20) ** 2, 1: 20 / 12, 2: (12 / 20) ** 2}
        for fld in src.image_fields:
            spec = next(iter(fld.spectra.values()))
            ratio = (spec(12000 * u.AA) / spec(20000 * u.AA)).value
            assert ratio == pytest.approx(expected[fld.meta["galaxy_id"]],
                                          rel=1e-3)
        src_all_late = _toy(n_cut=10.0)
        assert src_all_late.meta["n_early"] == 0

    def test_mag_limit_drops_faint_galaxies(self):
        assert _toy(mag_limit=24.0).meta["n_galaxies"] == 2

    def test_custom_catalogue_path(self):
        src = _toy()
        assert src.meta["catalogue"] == str(TOY_CATALOGUE)

    def test_empty_selection_raises(self):
        with pytest.raises(ValueError, match="zero galaxies"):
            _toy(mag_limit=10.0)

    def test_large_galaxy_is_coarsened_not_truncated(self, tmp_path):
        path = _write_single(tmp_path / "big.fits", r_eff=10.0, dx=0, dy=0)
        src = _toy(catalogue_path=path, pixel_scale=0.01, max_pixels=256)
        fld = src.image_fields[0]
        assert max(fld.data.shape) <= 257
        assert abs(fld.header["CDELT2"]) * 3600 > 0.01
        assert fld.data.sum() == pytest.approx(10 ** (-0.4 * 18.0), rel=1e-5)

    def test_clumps_conserve_total_flux(self):
        src_smooth = _vendored(fov=30.0)
        src_clumpy = _vendored(fov=30.0, clumps=True, seed=5)
        assert src_clumpy.meta["n_clumps"] > 0
        assert len(src_clumpy.table_fields) == 1
        assert _total_flux(src_clumpy) == pytest.approx(
            _total_flux(src_smooth), rel=1e-5)
        tbl = src_clumpy.table_fields[0].field
        assert set(tbl["ref"]).issubset(src_clumpy.table_fields[0].spectra)

    def test_seed_reproducible(self):
        a = _vendored(fov=30.0, clumps=True, seed=11)
        b = _vendored(fov=30.0, clumps=True, seed=11)
        c = _vendored(fov=30.0, clumps=True,
                      seed=np.random.default_rng(11))
        ta, tb, tc = (s.table_fields[0].field for s in (a, b, c))
        np.testing.assert_array_equal(ta["x"], tb["x"])
        np.testing.assert_array_equal(ta["weight"], tc["weight"])

    def test_global_rng_untouched(self):
        state = np.random.get_state()
        _vendored(fov=20.0, clumps=True, seed=1)
        _vendored(fov=20.0, clumps=True)
        after = np.random.get_state()
        np.testing.assert_array_equal(state[1], after[1])
        assert state[2:] == after[2:]

    def test_catalogue_file_not_modified(self, tmp_path):
        path = tmp_path / "toy.fits"
        shutil.copy(TOY_CATALOGUE, path)
        before = path.read_bytes()
        _toy(catalogue_path=path, clumps=True, seed=1)
        assert path.read_bytes() == before


# ---------------------------------------------------------------------------
# Real SEDs / filters and the optical train (network)
# ---------------------------------------------------------------------------

def _observe(src):
    from scopesim_templates.rc import load_example_optical_train
    opt = load_example_optical_train()
    opt["source_fits_keywords"].include = False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        opt.observe(src)
    return opt.image_planes[0].data, opt.cmds["!INST.pixel_scale"]


@pytest.mark.webtest
class TestHUDFWithRealData:
    def test_default_seds_magnitudes(self):
        """Image sum x spectrum reproduces the catalogue F160W mag."""
        src = hubble_ultra_deep_field(fov=10.0)
        fld = max(src.image_fields, key=lambda f: f.data.sum())
        spec = next(iter(fld.spectra.values()))
        from spextra import Passband
        mag0 = Observation(spec, Passband("HST/WFC3_IR.F160W")).effstim(
            u.ABmag).value
        mag = mag0 - 2.5 * np.log10(fld.data.sum(dtype=float))
        assert mag == pytest.approx(fld.meta["m_F160W"], abs=2e-3)

    def test_observed_field_has_signal_at_right_place(self):
        src = hubble_ultra_deep_field(fov=20.0, mag_limit=26.0)
        image, pixel_scale = _observe(src)
        median = np.median(image)
        assert image.max() > median * 1.005      # old bug: max == median
        brightest = max(src.image_fields, key=lambda f: f.data.sum())
        ny, nx = image.shape
        exp_x = (nx - 1) / 2 + brightest.header["CRVAL1"] * 3600 / pixel_scale
        exp_y = (ny - 1) / 2 + brightest.header["CRVAL2"] * 3600 / pixel_scale
        iy, ix = np.unravel_index(np.argmax(image), image.shape)
        assert abs(ix - exp_x) <= 1.5 and abs(iy - exp_y) <= 1.5

    def test_observed_orientation_and_flux(self, tmp_path):
        """PA = 45 deg (NE) -> major axis 45 deg CCW from +x (East) on the
        image plane; image flux equals that of point sources of the same
        magnitude."""
        path = _write_single(tmp_path / "one.fits", pa=45.0, dx=10, dy=5)
        src = hubble_ultra_deep_field(catalogue_path=path, ra_center=TOY_RA,
                                      dec_center=TOY_DEC)
        image, pixel_scale = _observe(src)
        image = image - np.median(image)
        cx, cy, angle = _moments(image)
        ny, nx = image.shape
        assert cx - (nx - 1) / 2 == pytest.approx(10 / pixel_scale, abs=1)
        assert cy - (ny - 1) / 2 == pytest.approx(5 / pixel_scale, abs=1)
        assert _axis_diff(angle, 45.0) < 5.0

        fld = src.image_fields[0]
        point = Source(x=[10.0], y=[5.0], ref=[0],
                       weight=[fld.data.sum(dtype=float)],
                       spectra=[next(iter(fld.spectra.values()))])
        image_pt, _ = _observe(point)
        image_pt = image_pt - np.median(image_pt)
        assert image.sum() == pytest.approx(image_pt.sum(), rel=0.02)
