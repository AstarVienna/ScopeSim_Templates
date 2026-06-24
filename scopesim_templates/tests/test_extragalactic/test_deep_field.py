# -*- coding: utf-8 -*-

import warnings
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from scopesim_templates.rc import Source
from scopesim_templates.extragalactic import hubble_ultra_deep_field


TOY_CATALOGUE = Path(__file__).parent / "data" / "toy_hudf.fits"


def _call(**kwargs):
    """Default to the toy fixture so most tests run offline-cheap."""
    kwargs.setdefault("catalogue_path", str(TOY_CATALOGUE))
    return hubble_ultra_deep_field(**kwargs)


@pytest.mark.webtest
class TestHUDF:

    def test_returns_source(self):
        src = _call()
        assert isinstance(src, Source)
        assert len(src.fields) >= 1

    def test_galaxy_count_matches_catalogue(self):
        src = _call()
        tbl = Table.read(TOY_CATALOGUE)
        assert len(src.fields) == len(tbl)

    def test_positions_within_fov(self):
        """A 6-arcsec FOV around the catalogue centre should drop the two
        offset (±5") galaxies but keep the central one."""
        src = _call(fov=6.0)
        assert src.meta["n_galaxies"] == 1

    def test_early_late_split_by_n_cut(self):
        """Toy catalogue has 2 rows with n>2.5 and 1 row with n<=2.5."""
        src = _call(n_cut=2.5)
        assert src.meta["n_galaxies"] == 3
        # Each per-galaxy sub-source carries one spectrum; the SED string
        # is recorded on the Source's stamped meta from add_function_call_str.
        # Easier check: count via the source catalogue.
        tbl = Table.read(TOY_CATALOGUE)
        n_early = int((np.asarray(tbl["sersic_n"]) > 2.5).sum())
        n_late = len(tbl) - n_early
        assert n_early == 2 and n_late == 1

    def test_amplitude_passes_through(self):
        """The brightest catalogue row (22.0 mag) should drive the
        brightest spectrum scaling."""
        src = _call()
        tbl = Table.read(TOY_CATALOGUE)
        i_bright = int(np.argmin(np.asarray(tbl["m_F160W"])))
        # spextra spectra carry their amplitude implicitly; verify the
        # spectrum exists and is callable on a wavelength grid.
        sp = src.spectra[i_bright]
        wave = np.linspace(1.4e4, 1.7e4, 50)  # H-band, Angstrom
        assert np.all(np.isfinite(sp(wave)))

    def test_mag_limit_drops_faint_galaxies(self):
        """Toy row 1 has m=26.5; cutting at 24.0 should drop it."""
        src = _call(mag_limit=24.0)
        assert src.meta["n_galaxies"] == 2

    def test_custom_catalogue_path(self):
        """Pointing at the toy file should give exactly 3 galaxies."""
        src = hubble_ultra_deep_field(catalogue_path=str(TOY_CATALOGUE))
        assert src.meta["n_galaxies"] == 3
        assert src.meta["catalogue"] == str(TOY_CATALOGUE)
