import time

import pytest
import numpy as np

from scopesim_templates.misc import nebula as neb
from scopesim_templates.misc.nebula import nebulosity


def _power_spectrum_slope(img, k_range=(0.01, 0.2)):
    """Fit the slope of the 2D power spectrum of a square image."""
    n = img.shape[0]
    power = np.abs(np.fft.fft2(img - img.mean())) ** 2
    ky, kx = np.meshgrid(np.fft.fftfreq(n), np.fft.fftfreq(n), indexing="ij")
    kk = np.hypot(kx, ky)
    mask = (kk > k_range[0]) & (kk < k_range[1])
    return np.polyfit(np.log(kk[mask]), np.log(power[mask]), 1)[0]


class TestNebulosity:
    @pytest.mark.parametrize("shape", [(256, 256), (300, 517), (17, 33),
                                       (1, 64)])
    @pytest.mark.parametrize("envelope", [None, "shell", "blob"])
    def test_returns_normalised_positive_array(self, shape, envelope):
        img = nebulosity(shape, envelope=envelope, seed=1)
        assert img.shape == shape
        assert img.dtype == np.float32
        assert np.all(np.isfinite(img))
        assert img.min() >= 0
        assert img.max() == pytest.approx(1)

    def test_dtype_can_be_set(self):
        img = nebulosity((64, 64), seed=1, dtype=np.float64)
        assert img.dtype == np.float64

    def test_same_seed_gives_same_image(self):
        img1 = nebulosity((128, 128), seed=42)
        img2 = nebulosity((128, 128), seed=42)
        assert np.array_equal(img1, img2)

    def test_different_seeds_give_different_images(self):
        img1 = nebulosity((128, 128), seed=1)
        img2 = nebulosity((128, 128), seed=2)
        assert not np.allclose(img1, img2)

    def test_accepts_generator_as_seed(self):
        img1 = nebulosity((128, 128), seed=np.random.default_rng(7))
        img2 = nebulosity((128, 128), seed=np.random.default_rng(7))
        assert np.array_equal(img1, img2)

    def test_default_is_a_new_realisation(self):
        assert not np.array_equal(nebulosity((64, 64)), nebulosity((64, 64)))

    def test_config_seed_makes_default_reproducible(self):
        from scopesim_templates.rc import __config__
        old = __config__["!random.seed"]
        try:
            __config__["!random.seed"] = 9001
            assert np.array_equal(nebulosity((64, 64)), nebulosity((64, 64)))
        finally:
            __config__["!random.seed"] = old

    def test_brightness_distribution_is_skewed(self):
        img = nebulosity((512, 512), seed=3).ravel()
        skew = ((img - img.mean()) ** 3).mean() / img.std() ** 3
        assert skew > 1

    def test_shell_envelope_is_brightest_on_the_ring(self):
        n = 512
        img = nebulosity((n, n), envelope="shell", seed=5,
                         envelope_kwargs={"ellipticity": 0, "roughness": 0,
                                          "asymmetry": 0})
        yy, xx = np.indices(img.shape)
        rr = np.hypot(yy - n / 2, xx - n / 2) / n
        ring = img[np.abs(rr - 0.35) < 0.03].mean()
        centre = img[rr < 0.1].mean()
        corners = img[rr > 0.6].mean()
        assert ring > 3 * centre
        assert ring > 10 * corners

    def test_array_envelope_is_applied(self):
        env = np.zeros((64, 64))
        env[:, 32:] = 1
        img = nebulosity((64, 64), envelope=env, seed=1)
        assert np.all(img[:, :32] == 0)
        assert np.all(img[:, 32:] > 0)

    def test_wrong_envelope_shape_raises(self):
        with pytest.raises(ValueError):
            nebulosity((64, 64), envelope=np.ones((32, 32)))

    def test_unknown_envelope_raises(self):
        with pytest.raises(ValueError):
            nebulosity((64, 64), envelope="donut")

    def test_bad_filament_strength_raises(self):
        with pytest.raises(ValueError):
            nebulosity((64, 64), filament_strength=1.5)

    def test_all_features_off_still_works(self):
        img = nebulosity((128, 128), filament_strength=0, warp_strength=0,
                         seed=1)
        assert np.all(np.isfinite(img))

    def test_2k_image_is_fast(self):
        start = time.perf_counter()
        nebulosity((2048, 2048), envelope="shell", seed=1)
        assert time.perf_counter() - start < 10


class TestPowerLawField:
    @pytest.mark.parametrize("slope", [2.0, 2.7, 3.5])
    def test_has_requested_spectral_slope(self, slope):
        n = 512
        shape = (n, n)
        white = neb._white_spectrum(shape, np.random.default_rng(0))
        field = neb._power_law_field(white, neb._kgrid(shape), shape, slope,
                                     1 / n, 0)
        assert _power_spectrum_slope(field) == pytest.approx(-slope, abs=0.3)

    def test_without_filaments_or_warp_keeps_slope_in_output(self):
        img = nebulosity((512, 512), filament_strength=0, warp_strength=0,
                         seed=0)
        assert _power_spectrum_slope(img) == pytest.approx(-2.7, abs=0.4)


class TestUpsample:
    def test_doubling_matches_periodic_linear_interpolation(self):
        from scipy import ndimage
        small = np.random.default_rng(0).random((64, 96), dtype=np.float32)
        big = neb._upsample(small, (128, 192))
        ref = ndimage.zoom(small, 2, order=1, mode="grid-wrap",
                           grid_mode=True)
        assert np.allclose(big, ref, atol=1e-5)

    @pytest.mark.parametrize("shape", [(256, 256), (301, 517), (33, 7)])
    def test_multiple_levels_give_right_shape_and_range(self, shape):
        small_shape = neb._downsampled_shape(shape, 4)
        small = np.random.default_rng(0).random(small_shape, dtype=np.float32)
        big = neb._upsample(small, shape)
        assert big.shape == shape
        assert big.min() >= small.min() - 1e-6
        assert big.max() <= small.max() + 1e-6

    def test_off_pyramid_shape_raises(self):
        with pytest.raises(ValueError):
            neb._upsample(np.ones((10, 10), dtype=np.float32), (64, 64))
