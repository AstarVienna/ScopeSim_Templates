# -*- coding: utf-8 -*-
"""
Off-tree generator for the vendored HUDF catalogue used by
``hubble_ultra_deep_field()``.

Pulls the 3D-HST GOODS-South master photometry (Skelton+ 2014,
CDS J/ApJS/214/24) and the van der Wel+ 2012 CANDELS structural
catalogue (J/ApJS/203/24), cross-matches by phot_id within the HUDF
footprint, applies quality + magnitude cuts, and writes a small FITS
table to ``data/hudf_catalogue.fits`` next to this script.

This script is **not** part of the public API. It runs once and the
resulting FITS is committed alongside the package. Re-run it only to
refresh the vendored catalogue.

Usage::

    python -m scopesim_templates.extragalactic._make_hudf_catalogue
"""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table


HUDF_CENTER = SkyCoord(ra=53.1625 * u.deg, dec=-27.7914 * u.deg, frame="icrs")
HUDF_RADIUS = 2.0 * u.arcmin            # cone radius covering the HUDF footprint
MAG_LIMIT = 28.0                        # m_F160W upper cap
OUT_PATH = Path(__file__).parent / "data" / "hudf_catalogue.fits"


def _fetch_skelton():
    """Fetch the 3D-HST cone around the HUDF centre (Skelton+ 2014)."""
    from astroquery.vizier import Vizier

    v = Vizier(columns=["**"], row_limit=-1, timeout=600)
    result = v.query_region(
        HUDF_CENTER, radius=HUDF_RADIUS, catalog="J/ApJS/214/24/3dhstall")
    if not result:
        raise RuntimeError("Empty result from Vizier for Skelton 3D-HST.")
    tbl = result[0]
    # The Skelton 3D-HST master catalog mixes all five CANDELS fields;
    # GOODS-S photometry has Field == 'goodss' (or similar). Filter by sky
    # alone is sufficient since the HUDF cone is uniquely in GOODS-S.
    return tbl


def _fetch_vdw_h_band():
    """Fetch the van der Wel+ 2012 F160W (H-band) structural fits."""
    from astroquery.vizier import Vizier

    v = Vizier(columns=["**"], row_limit=-1, timeout=600)
    result = v.query_region(
        HUDF_CENTER, radius=HUDF_RADIUS, catalog="J/ApJS/203/24")
    if not result:
        raise RuntimeError("Empty result from Vizier for van der Wel.")
    tbl = result[0]
    # Keep only F160W (H-band) rows.
    mask = np.asarray(tbl["F"]) == "H"
    return tbl[mask]


def _cross_match(skelton, vdw):
    """Spatial match (0.5" tolerance) — integer IDs differ between releases."""
    sky_sk = SkyCoord(skelton["RAJ2000"], skelton["DEJ2000"], unit="deg")
    sky_vd = SkyCoord(vdw["RAJ2000"], vdw["DEJ2000"], unit="deg")
    idx, sep, _ = sky_vd.match_to_catalog_sky(sky_sk)
    keep = sep < 0.5 * u.arcsec
    vd = vdw[keep]
    sk = skelton[idx[keep]]
    return sk, vd


def build_catalogue():
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    print(f"Fetching Skelton+2014 within {HUDF_RADIUS} of "
          f"{HUDF_CENTER.to_string('hmsdms')}...")
    skelton = _fetch_skelton()
    print(f"  {len(skelton)} rows")
    print("Fetching van der Wel+2012 (H-band) ...")
    vdw = _fetch_vdw_h_band()
    print(f"  {len(vdw)} rows")

    print("Cross-matching ...")
    sk, vd = _cross_match(skelton, vdw)
    print(f"  {len(sk)} matched rows")

    # Prefer Skelton's spec-z where available, fall back to peak photo-z.
    z = np.where(np.isfinite(np.asarray(sk["zsp"])) & (np.asarray(sk["zsp"]) > 0),
                 np.asarray(sk["zsp"]),
                 np.asarray(sk["zpk"]))
    z = np.where(np.isfinite(z) & (z > 0), z, 0.0)

    ra = np.asarray(vd["RAJ2000"], dtype=float)
    dec = np.asarray(vd["DEJ2000"], dtype=float)
    m_f160w = np.asarray(vd["mag"], dtype=float)
    r_eff = np.asarray(vd["r"], dtype=float)                 # arcsec
    sersic_n = np.asarray(vd["n"], dtype=float)
    q = np.asarray(vd["q"], dtype=float)
    pa = np.asarray(vd["PA"], dtype=float)                   # deg
    qflag = np.asarray(vd["Q"], dtype=int)

    mask = (qflag <= 1) & (m_f160w < MAG_LIMIT)
    mask &= np.isfinite(r_eff) & (r_eff > 0)
    mask &= np.isfinite(sersic_n) & (sersic_n > 0)
    mask &= np.isfinite(q) & (q > 0)

    print(f"  {mask.sum()} rows after quality + mag<{MAG_LIMIT} cuts")

    out = Table(
        data={
            "id": np.arange(int(mask.sum()), dtype=np.int32),
            "ra": ra[mask],
            "dec": dec[mask],
            "z": z[mask],
            "m_F160W": m_f160w[mask],
            "r_eff_arcsec": r_eff[mask],
            "sersic_n": sersic_n[mask],
            "axis_ratio": q[mask],
            "position_angle_deg": pa[mask],
        },
    )
    out["ra"].unit = u.deg
    out["dec"].unit = u.deg
    out["m_F160W"].unit = u.mag
    out["r_eff_arcsec"].unit = u.arcsec
    out["position_angle_deg"].unit = u.deg

    out.meta["origin"] = ("Skelton+2014 (J/ApJS/214/24) cross-matched with "
                          "van der Wel+2012 (J/ApJS/203/24).")
    out.meta["footprint"] = (
        f"cone r={HUDF_RADIUS.to_value(u.arcmin):.2f} arcmin around HUDF centre")
    out.meta["mag_limit"] = MAG_LIMIT
    out.meta["filter"] = "HST/WFC3.F160W"

    out.write(OUT_PATH, overwrite=True)
    print(f"Wrote {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.1f} kB)")
    return out


if __name__ == "__main__":
    build_catalogue()
