import numpy as np
from astropy import units as u
from scopesim_templates.basic import stars


class TestStars:
    def test_returns_the_spectra_scaled_to_ab_zero(self):
        spts = ["A0V", "G2V", "M9V", "G2V"]
        amplitudes = [1, 2, 3, 4] * u.mag
        x, y = np.arange(4), np.arange(4)
        filter_name = "Paranal/HAWKI.Ks"
        stars(spts, amplitudes, filter_name, x, y)

        print(stars)
