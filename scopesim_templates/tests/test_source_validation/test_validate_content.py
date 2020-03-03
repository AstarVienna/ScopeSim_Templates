import pytest

from astropy import units as u
from astropy.table import Table
from astropy.io.fits import ImageHDU, BinTableHDU
from synphot import SourceSpectrum, SpectralElement

from scopesim_templates.basic.stars import star
from scopesim_templates.basic.galaxy import spiral_two_component


### Add all initialied examples of sources to be tested to this list
SOURCE_LIST = [star(filter_name="Ks", amplitude=10*u.mag),
               star(filter_name="Paranal/HAWKI.J", amplitude=10*u.Jansky),
               spiral_two_component()
               ]


### Run all tests on each source indivdually
@pytest.mark.parametrize("src", SOURCE_LIST)
class TestContentOfFields:
    def test_is_fields_a_list(self, src):
        assert isinstance(src.fields, list)

    def test_fields_list_has_only_imagehdus_or_tables(self, src):
        for field in src.fields:
            assert isinstance(field, (Table, ImageHDU))

    def test_any_tables_have_correct_column_names(self, src):
        for field in src.fields:
            if isinstance(field, Table):
                req_colnames = ["x", "y", "ref", "weight"]
                assert all([col in field.colnames for col in req_colnames])

    def test_any_imagehdus_have_correct_header_keywords(self, src):
        for field in src.fields:
            if isinstance(field, ImageHDU):
                req_keys = ["SPEC_REF", "NAXIS", "NAXIS1", "NAXIS2",
                            "CUNIT1", "CUNIT2", "CTYPE1", "CTYPE2", "CDELT1",
                            "CDELT2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
                            ]
                assert all([key in field.header for key in req_keys])
