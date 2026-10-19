"""Regression tests for bugs found in the 2026-10 audit (docs/audit).

Each test pins one fixed bug; the docstring names the old behaviour.
"""

import warnings

import pytest
import numpy as np
from numpy import testing as npt
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from synphot import SourceSpectrum, ConstFlux1D

from scopesim_templates.rc import Source
from scopesim_templates.utils.general_utils import (function_call_str,
                                                    make_img_wcs_header)


FLAT_SED = SourceSpectrum(ConstFlux1D, amplitude=1)


def _peak_world_offset(hdu):
    """Return the (x, y) world offset in arcsec of the brightest pixel."""
    iy, ix = np.unravel_index(np.argmax(hdu.data), hdu.data.shape)
    wcs = WCS(hdu.header)
    world = wcs.wcs_pix2world([[ix, iy]], 0)[0]
    return world * 3600


class TestFunctionCallStr:
    def test_does_not_raise_on_too_many_positional_args(self):
        """Used to raise StopIteration from inside the decorator."""
        def func(a, b=1):
            pass
        call = function_call_str(func, (1, 2, 3))
        assert "func(" in call

    def test_expands_var_keyword(self):
        """**kwargs used to be reported as a single 'kwargs=...' entry."""
        def func(a, **kwargs):
            pass
        call = function_call_str(func, (1,), {"n": 4})
        assert call.endswith("func(a=1, n=4)")

    def test_includes_defaults(self):
        def func(a, b=2):
            pass
        assert function_call_str(func, (1,)).endswith("func(a=1, b=2)")


class TestMicado:
    def test_star_import_works(self):
        """micado.__all__ held function objects, so `import *` raised."""
        namespace = {}
        exec("from scopesim_templates.micado import *", namespace)
        assert "flatlamp" in namespace

    def test_pinhole_mask_wrapper_forwards_arguments(self):
        """micado.pinhole_mask() took no arguments and always raised."""
        from scopesim_templates import micado
        src = micado.pinhole_mask([0, 1], [0, 1], np.arange(1, 2, 0.01))
        assert isinstance(src, Source)
        assert src.meta["object"] == "pinhole_mask"

    def test_pinhole_mask_does_not_modify_input_waves(self):
        from scopesim_templates.micado.pinhole_masks import pinhole_mask
        waves = np.arange(1, 2, 0.01)
        pinhole_mask([0], [0], waves)
        assert not isinstance(waves, u.Quantity)

    def test_flatlamp_is_centred(self):
        """The centre was computed in arcsec, not pixels (off by 1/2 px)."""
        from scopesim_templates.micado.flatlamp import flatlamp
        hdu = flatlamp().fields[0].field
        ny, nx = hdu.data.shape
        assert hdu.header["CRPIX1"] == (nx + 1) / 2
        assert hdu.header["CRPIX2"] == (ny + 1) / 2
        # symmetric illumination pattern
        npt.assert_allclose(hdu.data, hdu.data[::-1, ::-1])

    def test_line_list_image_matches_wcs(self):
        """The image array was transposed relative to its WCS."""
        from scopesim_templates.micado import spectral_calibrations as sc
        hdu = sc.line_list(width=15., height=1.).fields[0].field
        ny, nx = hdu.data.shape
        assert nx > ny     # 15" wide, 1" high slit
        assert hdu.header["CRPIX1"] == pytest.approx((nx + 1) / 2)
        assert hdu.header["CRPIX2"] == pytest.approx((ny + 1) / 2)


class TestWcsAxisOrder:
    def test_make_img_wcs_header_takes_nx_ny(self):
        hdr = make_img_wcs_header(0, 0, 1, image_size=(30, 10))
        assert (hdr["CRPIX1"], hdr["CRPIX2"]) == (15.5, 5.5)

    def test_source_from_array_non_square_is_centred(self):
        """CRPIX1/2 were swapped for non-square arrays."""
        from scopesim_templates.misc.misc import source_from_array
        arr = np.zeros((11, 31))
        arr[5, 15] = 1          # the central pixel
        src = source_from_array(arr, FLAT_SED, pixel_scale=1.,
                                amplitude=None, filter_curve="V")
        npt.assert_allclose(_peak_world_offset(src.fields[0].field), 0,
                            atol=1e-6)


class TestNoInputMutation:
    @staticmethod
    def _hdu():
        hdu = fits.ImageHDU(data=np.ones((5, 5)))
        hdu.header.update({"CRVAL1": 12.3, "CRVAL2": 4.5, "CRPIX1": 3,
                           "CRPIX2": 3, "CDELT1": 1/3600, "CDELT2": 1/3600,
                           "CUNIT1": "deg", "CUNIT2": "deg",
                           "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
                           "BUNIT": "Jy"})
        return hdu

    def test_source_from_imagehdu_with_flux(self):
        from scopesim_templates.misc.misc import source_from_imagehdu_with_flux
        hdu = self._hdu()
        header_before = hdu.header.copy()
        source_from_imagehdu_with_flux(image_hdu=hdu, pixel_scale=0.1)
        assert hdu.header == header_before

    def test_source_from_imagehdu(self):
        from scopesim_templates.misc.misc import source_from_imagehdu
        hdu = self._hdu()
        source_from_imagehdu(hdu, filter_name="J")
        assert "SPEC_REF" not in hdu.header

    def test_source_from_file_closes_file_and_keeps_data(self, tmp_path):
        from scopesim_templates.misc.misc import source_from_file
        path = tmp_path / "img.fits"
        data = np.arange(25.).reshape(5, 5)
        fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data)]).writeto(path)
        source_from_file(path, pixel_scale=1., sed=FLAT_SED, amplitude=None,
                         filter_curve="V", cut=10)
        path.unlink()       # fails on Windows if the file is still open
        npt.assert_array_equal(data, np.arange(25.).reshape(5, 5))


class TestSpiralTwoComponent:
    def test_offset_is_applied(self):
        from scopesim_templates.extragalactic.galaxies import \
            spiral_two_component
        src = spiral_two_component(offset=(3, -2))
        hdr = src.fields[0].header
        assert hdr["CRVAL1"] * 3600 == pytest.approx(3)
        assert hdr["CRVAL2"] * 3600 == pytest.approx(-2)

    def test_fluxes_must_be_a_pair(self):
        from scopesim_templates.extragalactic.galaxies import \
            spiral_two_component
        with pytest.raises(ValueError):
            spiral_two_component(fluxes=(10, 12, 14))

    @pytest.mark.webtest
    @pytest.mark.parametrize("fluxes, system", [
        ((12, 14), "AB"),                       # floats are AB mags
        ((12 * u.ABmag, 14 * u.ABmag), "AB"),
        ((13 * u.mag, 11 * u.mag), "vega"),
    ])
    def test_fluxes_scale_each_component(self, fluxes, system):
        """`fluxes` used to be silently ignored."""
        from spextra import Spextrum
        from scopesim_templates.extragalactic.galaxies import \
            spiral_two_component
        src = spiral_two_component(fluxes=fluxes, filter_curve="Ks")
        for field, flux in zip(src.fields, fluxes):
            assert field.data.sum() == pytest.approx(1)
            spec = src.spectra[field.header["SPEC_REF"]]
            mag = Spextrum(modelclass=spec).get_magnitude(
                "Ks", system_name=system)
            assert mag.value == pytest.approx(getattr(flux, "value", flux),
                                              abs=0.01)

    @pytest.mark.webtest
    def test_jansky_fluxes(self):
        from synphot import Observation, units
        from spextra import Passband
        from scopesim_templates.extragalactic.galaxies import \
            spiral_two_component
        src_1 = spiral_two_component(fluxes=(1, 1) * u.Jy, filter_curve="V")
        src_2 = spiral_two_component(fluxes=(2, 1) * u.Jy, filter_curve="V")
        filt = Passband("V")
        ref = [src.spectra[src.fields[0].header["SPEC_REF"]]
               for src in (src_1, src_2)]
        ratio = (Observation(ref[1], filt).effstim(units.PHOTLAM)
                 / Observation(ref[0], filt).effstim(units.PHOTLAM))
        assert ratio.value == pytest.approx(2, rel=1e-3)


@pytest.mark.webtest
class TestStars:
    def test_does_not_modify_callers_spec_types(self):
        from scopesim_templates.stellar import stars
        spts = ["A0V"]
        stars("V", [10, 11] * u.mag, spts, [0, 1], [0, 1])
        assert spts == ["A0V"]

    def test_non_pyckles_library_weights_follow_amplitudes(self):
        """Stars of one type all pointed at the 1st star's scaled spectrum."""
        from scopesim_templates.stellar import stars
        src = stars("V", [10, 15] * u.mag, ["A0V", "A0V"], [0, 1], [0, 1],
                    library="kurucz")
        tbl = src.fields[0].field
        assert len(src.spectra) == 1
        npt.assert_allclose(tbl["weight"][1] / tbl["weight"][0], 10**-2)


@pytest.mark.webtest
class TestVikingFields:
    def test_load_stars_source_with_explicit_centre(self):
        """Passing ra/dec raised UnboundLocalError."""
        from scopesim_templates.micado.viking_fields import viking_fields as vf
        src = vf.load_stars_source(ra=0, dec=0)
        assert len(src.fields[0].field) > 1

    def test_galaxies_are_scaled_individually(self):
        """Every galaxy used the first galaxy's flux; ns= raised."""
        from astropy.io import ascii
        from scopesim_templates.micado.viking_fields import viking_fields as vf
        tbl = ascii.read(vf.DATA_DIR / "faintgal_1.dat")
        n = len(tbl)
        src = vf.load_galaxies_source("1", ns=np.full(n, 2.), seed=1,
                                      pixel_scale=0.02)
        flux = np.array([f.field.data.sum() for f in src.fields[:20]])
        expected = 10**(-0.4 * tbl["H_mag"][:20].astype(float))
        npt.assert_allclose(flux / flux[0], expected / expected[0],
                            rtol=0.05)

    def test_random_seed_alias_is_deprecated(self):
        from astropy.utils.exceptions import AstropyDeprecationWarning
        from scopesim_templates.micado.viking_fields import viking_fields as vf
        with pytest.warns(AstropyDeprecationWarning):
            vf.load_galaxies_source("1", random_seed=1, pixel_scale=0.04)


@pytest.mark.webtest
class TestGalaxySed:
    def test_galaxy_accepts_sourcespectrum(self):
        """A plain synphot SourceSpectrum raised AttributeError."""
        from scopesim_templates.extragalactic import galaxy
        from synphot import Empirical1D
        sed = SourceSpectrum(Empirical1D, points=[3000, 6000, 9000, 12000],
                             lookup_table=[1, 1, 1, 1])
        src = galaxy(sed=sed, z=0.5, amplitude=15, filter_curve="V")
        assert isinstance(src, Source)

    def test_galaxy3d_scales_spextrum_input(self):
        """z and amplitude were silently ignored for Spextrum input."""
        from spextra import Spextrum
        from scopesim_templates.extragalactic import galaxy3d
        sed = Spextrum("kc96/s0")
        src = galaxy3d(sed=sed, amplitude=12 * u.ABmag, filter_curve="V",
                       r_eff=1, pixel_scale=0.2, ngrid=2)
        spec = next(iter(src.spectra.values()))
        mag = Spextrum(modelclass=spec).get_magnitude("V", system_name="AB")
        # sub-sources are smoothed/shifted copies; the scale must be ~12 mag
        assert abs(mag.value - 12) < 0.5


class TestDeprecations:
    def test_micado_empty_sky_is_deprecated_but_works(self):
        from astropy.utils.exceptions import AstropyDeprecationWarning
        from scopesim_templates import micado
        with pytest.warns(AstropyDeprecationWarning, match="darkness"):
            src = micado.empty_sky()
        assert isinstance(src, Source)

    def test_source_from_image_is_deprecated_but_works(self):
        from scopesim_templates.misc.misc import source_from_image
        with pytest.warns(DeprecationWarning, match="source_from_array"):
            src = source_from_image(arr=np.ones((3, 3)), sed=FLAT_SED,
                                    pixel_scale=1., amplitude=None,
                                    filter_curve="V")
        assert isinstance(src, Source)


@pytest.mark.webtest
def test_star_records_its_own_call():
    from scopesim_templates.stellar import star
    src = star("V", 10)
    assert src.meta["object"] == "star"
    assert ".star(" in src.meta["function_call"]


class TestGalaxyPositionAngle:
    @pytest.mark.parametrize("theta", [0, 45, 90])
    def test_theta_is_in_degrees(self, theta):
        """Sersic2D got theta in radians while the velocity field used deg."""
        from scopesim_templates.extragalactic.galaxies import _galaxy_setup
        gal = _galaxy_setup(1 * u.arcsec, 10, 3, n=1, ellip=0.6, theta=theta)
        img = gal.intensity
        y, x = np.indices(img.shape)
        w = img / img.sum()
        dx, dy = x - (w * x).sum(), y - (w * y).sum()
        # orientation of the major axis from the second moments
        angle = 0.5 * np.degrees(np.arctan2(2 * (w * dx * dy).sum(),
                                            (w * dx**2).sum()
                                            - (w * dy**2).sum()))
        assert angle % 180 == pytest.approx(theta % 180, abs=1)
