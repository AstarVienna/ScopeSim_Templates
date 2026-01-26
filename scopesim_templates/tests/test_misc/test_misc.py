from os.path import exists
import pytest
import synphot
from pytest import raises
from astropy import units as u
from astropy.io import fits
import numpy as np
from matplotlib.colors import LogNorm
from scopesim.utils import figure_factory

from scopesim_templates.misc import misc
from scopesim_templates.tests.pyobjects import source_objects as so
from scopesim_templates.rc import Source, load_example_optical_train

METIS_FILTER_PATH = r"F:/Work/irdb/METIS/filters/TC_filter_H2O-ice.dat"

PLOTS = False


@pytest.fixture(scope="module")
def starting_star_imagehdu():
    """Create the same star as in starting.ipynb"""
    n = 100
    sigma = 5
    x, y = np.meshgrid(np.arange(n), np.arange(n))
    img = np.exp(-1 * (((x - n / 2) / sigma) ** 2 + ((y - n / 2) / sigma) ** 2))

    # Fits headers of the image. Yes it needs a WCS
    hdr = fits.Header({
        "NAXIS": 2,
        "NAXIS1": n,
        "NAXIS2": n,
        "CRPIX1": (n + 1) / 2,
        "CRPIX2": (n + 1) / 2,
        "CRVAL1": 0,
        "CRVAL2": 0,
        "CDELT1": 0.2 / 3600,
        "CDELT2": 0.2 / 3600,
        "CUNIT1": "DEG",
        "CUNIT2": "DEG",
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
    })

    # Creating an ImageHDU object
    hdu = fits.ImageHDU(data=img, header=hdr)
    return hdu


class TestSourceFromImageHDU:
    def test_initialises_for_basic_imagehdu(self):
        hdu = so._basic_imagehdu()
        src = misc.source_from_imagehdu(image_hdu=hdu, filter_name="J",
                                        pixel_unit_amplitude=20*u.mag)

        assert isinstance(src, Source)

    def test_raises_error_with_no_units_provided(self):
        hdu = so._basic_imagehdu()
        with raises(ValueError):
            misc.source_from_imagehdu(image_hdu=hdu, filter_name="J",
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

    # TODO: This test would be useful to check if supplying a filtername as a
    #       path does actually work, which is currently not tested anywhere
    #       else and potentially broken. Find a way to include a test file for
    #       this and then reenable the test!
    @pytest.mark.skipif(not exists(METIS_FILTER_PATH),
                        reason="Test only works on Kieran's local machine")
    def test_is_happy_with_metis_filter_and_bunit(self):
        hdu = so._basic_imagehdu()
        hdu.header["BUNIT"] = "Jy"
        src = misc.source_from_imagehdu(image_hdu=hdu,
                                        filter_name=METIS_FILTER_PATH)

        assert isinstance(src, Source)

    def test_actually_produces_image(self, starting_star_imagehdu):
        src1 = misc.source_from_imagehdu(
            starting_star_imagehdu,
            filter_name="J",
            pixel_unit_amplitude=10e12 * u.Jy)

        width_height = 4096
        opt = load_example_optical_train()
        opt['psf'].include = False
        opt.cmds["!OBS.psf_fwhm"] = 0.01
        opt.cmds["!TEL.area"] = 1000 * u.m**2
        opt.cmds["!INST.pixel_scale"] = 0.004
        opt.cmds["!INST.plate_scale"] = 0.4
        opt.cmds["!DET.width"] = width_height
        opt.cmds["!DET.height"] = width_height
        opt.cmds["!DET.dit"] = 30
        opt.cmds["!DET.ndit"] = 120

        opt["source_fits_keywords"].include = False

        opt.observe(src1)
        hdul = opt.readout()[0]

        data = hdul[1].data
        # Is the background okay?
        assert 500 < np.median(data) < 5000
        # Is there a bright source?
        assert data.max() > 10000000
        # Is the bright source approximately in the center.
        x_cen, y_cen = np.unravel_index(data.argmax(), data.shape)
        # TODO: Figure out why source is not in the center!
        fudge = 100
        assert width_height / 2 - fudge < x_cen < width_height / 2 + fudge
        assert width_height / 2 - fudge < y_cen < width_height / 2 + fudge

        if PLOTS:
            fig, ax = figure_factory(3, 1)
            ax[0].imshow(src1.fields[0].data)
            ax[1].imshow(opt.image_planes[0].data, norm=LogNorm())
            ax[2].imshow(hdul[1].data, norm=LogNorm())
            fig.show()

    def test_scales(self, starting_star_imagehdu):
        """Test whether source_from_imagehdu scales w/ pixel_unit_amplitude."""
        flux1 = 10e12 * u.Jy
        src1 = misc.source_from_imagehdu(
            starting_star_imagehdu,
            filter_name="J",
            pixel_unit_amplitude=flux1)

        src2 = misc.source_from_imagehdu(
            starting_star_imagehdu,
            filter_name="J",
            pixel_unit_amplitude=flux1 * 10)

        width_height = 4096
        opt = load_example_optical_train()
        opt['psf'].include = False
        opt.cmds["!OBS.psf_fwhm"] = 0.01
        opt.cmds["!TEL.area"] = 1000 * u.m**2
        opt.cmds["!INST.pixel_scale"] = 0.004
        opt.cmds["!INST.plate_scale"] = 0.4
        opt.cmds["!DET.width"] = width_height
        opt.cmds["!DET.height"] = width_height
        opt.cmds["!DET.dit"] = 30
        opt.cmds["!DET.ndit"] = 120

        opt["source_fits_keywords"].include = False

        opt.observe(src1)
        hdul1 = opt.readout()[0]
        data1 = hdul1[1].data

        opt.observe(src2, update=True)
        hdul2 = opt.readout()[0]
        data2 = hdul2[1].data

        if PLOTS:
            fig, ax = figure_factory(1, 1)
            ax.imshow(data2 / data1)
            fig.show()

        max1 = data1.max()
        max2 = data2.max()
        assert 9 < max2 / max1 < 11


class TestPointSource:
    @pytest.mark.webtest
    def test_initialize_from_spextra(self):
        src = misc.point_source("pickles/a0v", amplitude=16)

        assert isinstance(src, Source)

    def test_initialize_from_synphot(self):
        sp = synphot.SourceSpectrum(synphot.Empirical1D, points=[1000] * 4,
                                    lookup_table=[1] * 4)
        src = misc.point_source(sed=sp)
        assert isinstance(src, Source)


class TestUniformSource:
    @pytest.mark.webtest
    def test_initialize_from_spextra(self):
        src = misc.uniform_source("pickles/a0v", amplitude=16)

        assert isinstance(src, Source)

    def test_initialize_from_synphot(self):
        sp = synphot.SourceSpectrum(synphot.Empirical1D, points=[1000] * 4,
                                    lookup_table=[1] * 4)
        src = misc.uniform_source(sed=sp)
        assert isinstance(src, Source)


def test_poorman_cube_source_is_working():
    cube = so._make_dummy_cube(scale=0.2, wave_unit=u.AA, ref_wave=5000,
                               wave_step=1, wave_type="WAVE",
                               bunit="erg / (s cm2 Angstrom)")

    hdul = fits.HDUList(cube)
    cube_source = misc.poorman_cube_source(hdu=hdul, ext=0)

    assert isinstance(cube_source, Source)
