from pathlib import Path

import numpy as np
from astropy.io import ascii
from astropy import units as u
import yaml

import pyckles

from ...stellar import stars
from ...extragalactic import elliptical
from ...rc import ter_curve_utils as tcu
from ...rc import Source


DATA_DIR = Path(__file__).parent


def viking_field(star_cat_id="illum", gal_cat_id="1", pixel_scale=0.004,
                 ra=None, dec=None, gal_field_size=60, gal_x=None, gal_y=None,
                 random_seed=None):
    """
    Make a field with background stars and galaxies according to Ric's input.

    Any combination of the star and galaxy catalogues is possible.


    Parameters
    ----------
    star_cat_id : str
        Default "illum". Options ["illum", "science", "stdstar"]

    gal_cat_id : str
        Default "1". Options ["1", "2", "3"]

    pixel_scale : float
        [arcsec] Default 0.004"

    ra, dec : float, None
        [deg] Default None. Coords for centre or star field.
        If None, field centre is taken as average of all star positions

    gal_field_size : float
        [arcsec] Default 60". Size of galaxy field for randomly assigning coords

    gal_x, gal_y : list of floats
        [arcsec] Default None.
        If None, galaxy coordinates are randomly assigned within a 1'x1' field

    random_seed : int, None
        For reproducible galaxy positions


    Returns
    -------
    viking_field : Source
        Combined source of all stars and galaxies


    Examples
    --------
    ::

        from scopesim_templates.calibration.micado import viking_fields as vf

        # Standard random star and galaxy generation at 0.004" resolution
        src = vf.viking_field(star_cat_id="illum", gal_cat_id="1")

        # For oversampling the galaxies and fixing their positions:
        src = vf.viking_field(star_cat_id="science", gal_cat_id="2",
                              pixel_scale=0.002, random_seed=9001)

    """
    star_src = load_stars_source(cat_id=star_cat_id, ra=ra, dec=dec)
    gal_src = load_galaxies_source(cat_id=gal_cat_id, pixel_scale=pixel_scale,
                                   box_size=gal_field_size, xs=gal_x, ys=gal_y)

    return gal_src + star_src


def load_stars_source(cat_id="illum", ra=None, dec=None):
    """
    Make a Source for one of Ric's Viking star fields.

    Parameters
    ----------
    cat_id : str
        Default "illum". Options ["illum", "science", "stdstar"]
    ra, dec : float, None
        [deg] Coords for centre or star field.
        If None, (RA, Dec) centre is taken from file header

    Returns
    -------
    src : Source

    """
    tbl = ascii.read(DATA_DIR / f"cat_{cat_id}.dat")
    tbl_meta = yaml.full_load("\n".join(tbl.meta["comments"]))

    if ra is None:
        ra0 = tbl_meta["RA_centre"]
    if dec is None:
        dec0 = tbl_meta["Dec_centre"]
    x, y = 3600 * (tbl["RA"] - ra0), 3600 * (tbl["Dec"] - dec0)

    mags = tbl["Hmag"].data
    src = stars(filter_name="H", amplitudes=mags, spec_types="A0V", x=x, y=y)

    return src


def load_galaxies_source(cat_id="1", pixel_scale=0.004, xs=None, ys=None,
                         angles=None, ellipticitys=None, ns=None, box_size=60,
                         random_seed=None):
    """
    Make a Source full of faint galaxies for one of Ric's background fields.

    .. warning: This funciton uses astropy's Sersic2D class, which can
       have unexpected consequences with GalFit (speak to Carmelo Archidiacono)


    Parameters
    ----------
    cat_id : str
        Default "1". Options ["1", "2", "3"]

    pixel_scale : float
        [arcsec] Default 0.004"

    xs, ys : list of floats
        [arcsec] : central positions of galaxies relative to field centre

    angles : list of floats
        [deg] rotation of galaxies

    ellipticitys : list of floats
        [0..1] 0=circle, 1=line

    ns : list of floats
        Sersic indicies

    box_size : float
        [arcsec] Default 60. Size of box for randomly assigned x,y positions

    random_seed : int, None
        For reproducible galaxy positions


    Returns
    -------
    gals_src : Source


    """
    if random_seed is not None:
        np.random.seed(random_seed)

    tbl = ascii.read(DATA_DIR / f"faintgal_{cat_id}.dat")

    h_mag = tbl["H_mag"].data.astype(float)             # [Vega mag]
    reff = tbl["Reff"].data.astype(float)               # [mas]
    gal_type = tbl["type"].data.astype(float)           # [0=ellip, 1=spiral]

    scale_factors = 10**(-0.4 * h_mag)
    r_effs = 1e-3 * reff                                # [arcsec]
    widths = 2 * r_effs

    # Or should there be a range > involved? E.g. Spiral 1-2, Ellip 3-5?
    # Ellip type 0, spiral type 1
    # Do spirals first
    if ns is None:
        n_gals = len(gal_type)
        ns = np.random.random(n_gals)   # Sersic index
        ns[gal_type == 1] += 1            # Spirals n=[1..2]
        ns[gal_type == 0] *= 2            # Ellips n=[3..5] => ([0..1] * 2) + 3
        ns[gal_type == 0] += 3

    if angles is None:
        angles = np.random.random(n_gals) * 360         # [deg]

    if ellipticitys is None:
        e_max = 0.85                                    # 0 = circular
        ellipticitys = e_max * np.random.random(n_gals)
        ellipticitys[gal_type == 0] = 0

    if xs is None:
        xs = np.random.random(n_gals) * box_size - 0.5 * box_size   # [arcsec]

    if ys is None:
        ys = np.random.random(n_gals) * box_size - 0.5 * box_size   # [arcsec]

    brown_lib = pyckles.SpectralLibrary("brown", return_style="synphot")
    spectrum = brown_lib["NGC_0584"]
    spectrum = tcu.scale_spectrum(spectrum, filter_name="H", amplitude=0*u.mag)

    srcs = []
    for i in range(n_gals):
        src = elliptical(r_eff=r_effs[i],
                         pixel_scale=pixel_scale,
                         filter_name="H",
                         amplitude=1,
                         spectrum=spectrum,
                         rescale_spectrum=False,
                         n=ns[i],
                         angle=angles[i],
                         ellipticity=ellipticitys[i],
                         width=widths[i],
                         height=widths[i],
                         normalization="total",
                         )
        src.fields[0].data *= scale_factors[0]
        src.shift(dx=xs[i], dy=ys[i])
        srcs += [src]

    gals_src = Source()
    gals_src.fields = [src.fields[0] for src in srcs]
    gals_src.spectra = [spectrum]
    gals_src._meta_dicts = [src.meta for src in srcs]

    return gals_src
