# -*- coding: utf-8 -*-
"""Hubble Ultra Deep Field analogue, built from a vendored HUDF catalogue."""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from ..rc import Source
from ..utils.general_utils import add_function_call_str
from .galaxies import galaxy


__all__ = ["hubble_ultra_deep_field"]


_DEFAULT_CATALOGUE = Path(__file__).parent / "data" / "hudf_catalogue.fits"

HUDF_RA = 53.1625    # deg
HUDF_DEC = -27.7914  # deg


def _select_rows(tbl, ra_center, dec_center, fov, mag_limit, mag_col):
    """Return a boolean mask for the square FOV + magnitude cut."""
    mask = np.ones(len(tbl), dtype=bool)
    if mag_limit is not None:
        mask &= np.asarray(tbl[mag_col]) <= float(mag_limit)
    if fov is not None:
        half = 0.5 * float(fov) * u.arcsec
        ctr = SkyCoord(ra=ra_center * u.deg, dec=dec_center * u.deg)
        rows = SkyCoord(ra=np.asarray(tbl["ra"]) * u.deg,
                        dec=np.asarray(tbl["dec"]) * u.deg)
        dra, ddec = ctr.spherical_offsets_to(rows)
        mask &= (np.abs(dra) <= half) & (np.abs(ddec) <= half)
    return mask


@add_function_call_str
def hubble_ultra_deep_field(pixel_scale=0.06,
                            fov=None,
                            ra_center=HUDF_RA,
                            dec_center=HUDF_DEC,
                            mag_limit=None,
                            filter_name="HST/WFC3_IR.F160W",
                            sed_early="brown/NGC0584",
                            sed_late="brown/NGC4254",
                            n_cut=2.5,
                            catalogue_path=None):
    """
    Source of an HUDF-analogue galaxy field.

    Loads a vendored catalogue of HUDF galaxies (Skelton+2014 photometry
    cross-matched with van der Wel+2012 F160W structural fits), trims it
    to a square FOV and magnitude limit, and assembles one
    ``scopesim.Source`` by summing the per-galaxy Sersic profiles produced
    by :func:`scopesim_templates.extragalactic.galaxy`.

    Each galaxy gets one of two user-supplied spextra SEDs based on a
    Sersic-index cut (``n > n_cut`` is treated as early-type), redshifted
    to the catalogue ``z`` and scaled to the catalogue ``m_F160W`` in
    ``filter_name``.

    Parameters
    ----------
    pixel_scale : float
        Output pixel scale in arcsec/pixel. Default 0.06 (HST WFC3 NIR).
    fov : float, optional
        Side length of the square FOV in arcsec. ``None`` (default) uses
        the full extent of the vendored catalogue (~4 arcmin diameter
        around the HUDF centre).
    ra_center, dec_center : float
        FOV centre, degrees. Default is the HUDF pointing.
    mag_limit : float, optional
        Drop galaxies fainter than this F160W AB mag. ``None`` uses the
        catalogue's own limit (28.0 mag for the default file).
    filter_name : str
        SVO filter ID used to scale the SEDs. Default ``"HST/WFC3_IR.F160W"``
        matches the GALFIT magnitude column in the vendored catalogue.
    sed_early, sed_late : str
        spextra SED names for ``n > n_cut`` and ``n <= n_cut`` galaxies.
        Defaults are Brown+2014 archetypes: ``NGC0584`` (E0) and
        ``NGC4254`` (Sc spiral). The kc96 library is *not* suitable
        here because it only extends to ~1 um.
    n_cut : float
        Sersic-index threshold separating early- and late-type SED
        assignment.
    catalogue_path : str or Path, optional
        Override the vendored FITS catalogue. The file must contain the
        same columns: ``ra, dec, z, m_F160W, r_eff_arcsec, sersic_n,
        axis_ratio, position_angle_deg``.

    Returns
    -------
    scopesim.Source

    Examples
    --------
    >>> from scopesim_templates.extragalactic import hubble_ultra_deep_field
    >>> src = hubble_ultra_deep_field(fov=20.0, mag_limit=25.0)
    """
    path = Path(catalogue_path) if catalogue_path is not None else _DEFAULT_CATALOGUE
    tbl = Table.read(path)

    mask = _select_rows(tbl, ra_center, dec_center, fov, mag_limit, "m_F160W")
    sel = tbl[mask]
    if len(sel) == 0:
        raise ValueError(
            f"hubble_ultra_deep_field: zero galaxies after fov={fov} and "
            f"mag_limit={mag_limit} cuts on {path}."
        )

    srcs = []
    for row in sel:
        n = float(row["sersic_n"])
        sed = sed_early if n > n_cut else sed_late
        q = float(row["axis_ratio"])
        srcs.append(galaxy(
            sed=sed,
            z=float(row["z"]),
            amplitude=float(row["m_F160W"]) * u.ABmag,
            filter_curve=filter_name,
            pixel_scale=pixel_scale,
            r_eff=float(row["r_eff_arcsec"]),
            n=n,
            ellip=1.0 - q,
            theta=float(row["position_angle_deg"]),
            ra=float(row["ra"]),
            dec=float(row["dec"]),
        ))

    src = sum(srcs[1:], start=srcs[0])
    src.meta.update({
        "object": "Hubble Ultra Deep Field analogue",
        "ra_center": ra_center,
        "dec_center": dec_center,
        "fov": fov,
        "mag_limit": mag_limit,
        "filter_name": filter_name,
        "sed_early": sed_early,
        "sed_late": sed_late,
        "n_cut": n_cut,
        "n_galaxies": len(sel),
        "catalogue": str(path),
    })
    return src
