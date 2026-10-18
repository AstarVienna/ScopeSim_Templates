"""Observe every public template with ScopeSim's example optical train.

Building a ``Source`` is not enough: the HUDF template once produced valid
Sources whose light never reached the detector (absolute RA/Dec in CRVAL).
These tests check that each template's flux actually lands on the image plane.
"""

import pytest
import numpy as np
from astropy import units as u
from synphot import SourceSpectrum, ConstFlux1D

import scopesim_templates as st
from scopesim_templates import stellar, extragalactic, misc, calibration
from scopesim_templates.rc import Source, load_example_optical_train


SED = SourceSpectrum(ConstFlux1D, amplitude=1)

BUILDERS = {
    "star": lambda: stellar.star("V", 12),
    "stars": lambda: stellar.stars("V", [12, 13] * u.mag, ["A0V", "G2V"],
                                   [-2, 2], [0, 1]),
    "star_field": lambda: stellar.star_field(10, 12, 15, 20, seed=1),
    "star_grid": lambda: stellar.star_grid(9, 12, 15, separation=2),
    "cluster": lambda: stellar.cluster(mass=200, distance=2000,
                                       core_radius=0.05, seed=1),
    "globular_cluster": lambda: stellar.globular_cluster(
        density=0.05, fov=20, distance_modulus=13.6, seed=1),
    "galactic_centre": lambda: stellar.galactic_centre("2024-01-01"),
    "s301": lambda: stellar.s301("J2031.8"),
    "omega_cen": lambda: stellar.omega_cen(fov=3, seed=1),
    # NB: kc96 SEDs stop at ~1 um, i.e. outside the example train's NIR band
    "galaxy": lambda: extragalactic.galaxy("brown/NGC0584", amplitude=12,
                                           filter_curve="V", r_eff=2,
                                           pixel_scale=0.2),
    "elliptical": lambda: extragalactic.elliptical(2, 0.2, "V", 12, width=20,
                                                   height=20),
    "spiral_two_component": lambda: extragalactic.spiral_two_component(
        extent=30 * u.arcsec),
    "hubble_ultra_deep_field": lambda: st.hubble_ultra_deep_field(
        fov=60),
    "point_source": lambda: misc.point_source(SED, amplitude=12),
    "nebulosity": lambda: misc.source_from_array(
        misc.nebulosity((64, 64), seed=1), SED, pixel_scale=0.3,
        amplitude=10, filter_curve="V"),
}


@pytest.fixture(scope="module")
def optical_train():
    opt = load_example_optical_train()
    opt["source_fits_keywords"].include = False
    return opt


@pytest.mark.webtest
@pytest.mark.parametrize("name", list(BUILDERS))
def test_template_flux_reaches_image_plane(name, optical_train):
    src = BUILDERS[name]()
    assert isinstance(src, Source)
    optical_train.observe(src, update=True)
    image = optical_train.image_planes[0].data
    assert np.all(np.isfinite(image))
    # Real signal: the image plane is noiseless, so a flat image (max ==
    # median) means no flux arrived. HUDF galaxies are faint vs the sky.
    assert image.max() > 1.01 * np.median(image), name
