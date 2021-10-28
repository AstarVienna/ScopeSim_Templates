from astropy import units as u
from astropy.io import fits
from scopesim import Source
from scopesim_templates.misc.misc import poorman_cube_source
from scopesim_templates.tests.pyobjects.source_objects import _make_dummy_cube

from scopesim_templates.rc import Source


class TestSourceFromImageHDU:
    def test_




def test_poorman_cube_source_is_working():
    cube = _make_dummy_cube(scale=0.2, wave_unit=u.AA, ref_wave=5000,
                            wave_step=1, wave_type="WAVE", bunit="erg / (s cm2 Angstrom)")

    hdul = fits.HDUList(cube)
    cube_source = poorman_cube_source(hdu=hdul, ext=0)

    isinstance(cube_source, Source)