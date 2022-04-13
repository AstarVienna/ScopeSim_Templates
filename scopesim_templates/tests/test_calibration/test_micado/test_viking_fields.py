import numpy as np
from astropy.io import fits
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scopesim_templates.calibration.micado import viking_fields as vf


class TestVikingCatalogues:
    def test_load_galaxy_cat(self):
        srcs = vf.load_galaxy_sources()

        plt.figure(figsize=(20,20))
        j = np.random.randint(0, 1131, 26)
        for i in range(1, 26):
            plt.subplot(5,5,i)
            plt.imshow(srcs[j[i]].fields[0].data, norm=LogNorm())

        plt.show()
        pass
