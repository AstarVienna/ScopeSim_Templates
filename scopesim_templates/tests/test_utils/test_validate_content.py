import pytest

import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.io.fits import ImageHDU
from synphot import SourceSpectrum

from scopesim_templates.stellar import star
from scopesim_templates.extragalactic import galaxy, galaxy3d, spiral_two_component
from scopesim_templates.misc.misc import source_from_array


@pytest.fixture(scope="module")
def source_list():
    """All initialized examples of sources to be tested"""
    return [
        star(filter_name="Ks", amplitude=10*u.mag),
        star(filter_name="Paranal/HAWKI.J", amplitude=10*u.Jansky),
        spiral_two_component(),
        galaxy(sed="kc96/s0"),
        galaxy3d(sed="kc96/s0", ngrid=10),
        source_from_array(
            arr=np.ones(shape=(100, 100)),
            sed="kc96/s0",
            amplitude=15,
            pixel_scale=0.2,
            filter_curve="g",
        ),
    ]


# Run all tests on each source indivdually
@pytest.mark.webtest
@pytest.mark.usefixtures("source_list")
class TestContentOfFields:
    """Test the <Source>.fields attribute to check that the content is correct."""

    def test_is_fields_a_list(self, source_list):
        for src in source_list:
            assert isinstance(src.fields, list)

    def test_fields_list_has_only_imagehdus_or_tables(self, source_list):
        for src in source_list:
            for field in src.fields:
                if not isinstance(field, Table):
                    assert isinstance(field.field, (Table, ImageHDU))

    def test_any_tables_have_correct_column_names(self, source_list):
        for src in source_list:
            for field in src.fields:
                if not isinstance(field, Table):
                    continue
                req_colnames = ["x", "y", "ref", "weight"]
                assert all(col in field.field.colnames for col in req_colnames)

    def test_any_imagehdus_have_correct_header_keywords(self, source_list):
        for src in source_list:
            for field in src.fields:
                if isinstance(field.field, ImageHDU):
                    req_keys = ["SPEC_REF", "NAXIS", "NAXIS1", "NAXIS2",
                                "CUNIT1", "CUNIT2", "CTYPE1", "CTYPE2",
                                "CDELT1", "CDELT2", "CRVAL1", "CRVAL2",
                                "CRPIX1", "CRPIX2"]
                    assert all(key in field.header for key in req_keys)


@pytest.mark.webtest
class TestContentOfSpectra:
    """Test the <Source>.spectra attribute to check that all spectra are useful."""

    def test_is_spectra_a_list_or_dict(self, source_list):
        for src in source_list:
            # ScopeSim pre-0.9 is a list, afterwards dict. Allow both here.
            assert isinstance(src.spectra, (list, dict))

    def test_spectra_list_has_only_synphot_sourcespectrum_objects(self, source_list):
        for src in source_list:
            for spectrum in src.spectra.values():
                assert isinstance(spectrum, SourceSpectrum)


@pytest.mark.webtest
class TestConnectionBetweenFieldsAndSpectra:
    """Test that all spectra referenced by .fields are in .spectra."""

    def test_all_spectra_in_table_ref_column_exist(self, source_list):
        for src in source_list:
            for field in src.fields:
                if isinstance(field.field, Table):
                    for ref in field.field["ref"]:
                        assert isinstance(src.spectra[ref], SourceSpectrum)

    def test_all_spectra_refereced_in_imagehdu_header_exist(self, source_list):
        for src in source_list:
            for field in src.fields:
                if isinstance(field.field, ImageHDU):
                    ref = field.header["SPEC_REF"]
                    assert isinstance(src.spectra[ref], SourceSpectrum)
