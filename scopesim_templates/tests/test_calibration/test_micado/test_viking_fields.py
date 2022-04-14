import pytest

import numpy as np
from astropy.io import fits
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scopesim_templates.calibration.micado import viking_fields as vf

PLOTS = False


class TestVikingCatalogues:
    @pytest.mark.parametrize("cat_id", [("1"), ("2"), ("3")])
    def test_load_galaxies_source(self, cat_id):
        gal_src = vf.load_galaxies_source(cat_id)

        if PLOTS:
            plt.figure(figsize=(20, 20))
            n = 10
            for i, r in enumerate(np.random.randint(0, 1131, n*n)):
                plt.subplot(n, n, i+1)
                plt.imshow(gal_src.fields[r].data, norm=LogNorm())
            plt.show()

        assert len(gal_src.fields) > 1
        assert len(gal_src.spectra) == 1
        assert len(gal_src.fields) == len(gal_src._meta_dicts)

    @pytest.mark.parametrize("cat_id",
                             [("stdstar"),
                              ("illum"),
                              ("science")])
    def test_load_stars_source(self, cat_id):
        star_src = vf.load_stars_source(cat_id)

        if PLOTS:
            plt.figure(figsize=(10, 10))
            x, y = star_src.fields[0]["x"], star_src.fields[0]["y"]
            size = 2.5*np.log10(star_src.fields[0]["weight"])
            size -= min(size)
            plt.scatter(x, y, s=size)
            plt.show()

        assert len(star_src.fields) == 1
        assert len(star_src.fields[0]) > 1
        assert len(star_src.spectra) == 1
        assert len(star_src.fields) == len(star_src._meta_dicts)

    def test_load_viking_field(self):
        viking_src = vf.viking_field()

        gals = viking_src.image_fields
        gal_x = [gal.header["CRVAL1"] * 3600 for gal in gals]
        gal_y = [gal.header["CRVAL2"] * 3600 for gal in gals]

        stars = viking_src.table_fields
        star_x, star_y = stars[0]["x"].data, stars[0]["y"].data

        if not PLOTS:
            plt.figure(figsize=(10, 10))
            plt.scatter(gal_x, gal_y, c="b")
            plt.scatter(star_x, star_y, c="r")
            plt.show()

        assert len(viking_src.fields) > 1
        assert len(viking_src.fields) == len(gals) + len(stars)
