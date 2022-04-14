from os.path import exists
import pytest
import synphot
from pytest import raises
from astropy import units as u
from astropy.io import fits

from scopesim_templates.misc import misc
from scopesim_templates.tests.pyobjects import source_objects as so
from scopesim_templates.rc import Source

METIS_FILTER_PATH = r"F:/Work/irdb/METIS/filters/TC_filter_H2O-ice.dat"


class TestSourceFromImageHDU:
    def test_initialises_for_basic_imagehdu(self):
        hdu = so._basic_imagehdu()
        src = misc.source_from_imagehdu(image_hdu=hdu, filter_name="J",
                                        pixel_unit_amplitude=20*u.mag)

        assert isinstance(src, Source)

    def test_raises_error_with_no_units_provided(self):
        hdu = so._basic_imagehdu()
        with raises(ValueError):
            src = misc.source_from_imagehdu(image_hdu=hdu, filter_name="J",
                                            pixel_unit_amplitude=20)

    def test_is_happy_with_only_bunit(self):
        hdu = so._basic_imagehdu()
        hdu.header["BUNIT"] = "erg s-1 cm-2 um-1"
        src = misc.source_from_imagehdu(image_hdu=hdu, filter_name="J")

        assert isinstance(src, Source)

    def test_is_happy_with_spanish_vo_filter_name(self):
        hdu = so._basic_imagehdu()
        filter_name = "Generic/Johnson_UBVRIJHKL.N"
        src = misc.source_from_imagehdu(image_hdu=hdu, filter_name=filter_name,
                                        pixel_unit_amplitude=20 * u.Jy)

        assert isinstance(src, Source)

    @pytest.mark.skipif(not exists(METIS_FILTER_PATH),
                        reason="Test only works on Kieran's local machine")
    def test_is_happy_with_metis_filter_and_bunit(self):
        hdu = so._basic_imagehdu()
        hdu.header["BUNIT"] = "Jy"
        src = misc.source_from_imagehdu(image_hdu=hdu,
                                        filter_name=METIS_FILTER_PATH)

        assert isinstance(src, Source)


class TestPointSource:
    def test_initialize_from_spextra(self):
        src = misc.point_source("pickles/a0v", amplitude=16)

        assert isinstance(src, Source)

    def test_initialize_from_synphot(self):
        sp = synphot.SourceSpectrum(synphot.Empirical1D, points=[1000, 10000], lookup_table=[1, 1])
        src = misc.point_source(sed=sp)
        assert isinstance(src, Source)


class TestUniformSource:
    def test_initialize_from_spextra(self):
        src = misc.uniform_source("pickles/a0v", amplitude=16)

        assert isinstance(src, Source)

    def test_initialize_from_synphot(self):
        sp = synphot.SourceSpectrum(synphot.Empirical1D, points=[1000, 10000], lookup_table=[1, 1])
        src = misc.uniform_source(sed=sp)
        assert isinstance(src, Source)


def test_poorman_cube_source_is_working():
    cube = so._make_dummy_cube(scale=0.2, wave_unit=u.AA, ref_wave=5000,
                            wave_step=1, wave_type="WAVE", bunit="erg / (s cm2 Angstrom)")

    hdul = fits.HDUList(cube)
    cube_source = misc.poorman_cube_source(hdu=hdul, ext=0)

    assert isinstance(cube_source, Source)

