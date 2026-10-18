# -*- coding: utf-8 -*-
"""
Off-tree generator for the vendored HUDF catalogue used by
``scopesim_templates.extragalactic.hubble_ultra_deep_field``.

Queries VizieR for

* **van der Wel et al. 2012, ApJS 203, 24** (``J/ApJS/203/24/table2``):
  GALFIT single-Sersic F160W fits (magnitude, half-light semi-major axis,
  Sersic index, axis ratio, position angle), GOODS-S, quality flag ``Q = 0``;
* **Skelton et al. 2014, ApJS 214, 24** (``J/ApJS/214/24/3dhstall``):
  3D-HST photometry (spectroscopic / peak photometric redshift, FAST stellar
  mass, F606W flux),

cross-matches them positionally (0.5 arcsec), and writes a small FITS table
to ``scopesim_templates/extragalactic/data/hudf_catalogue.fits``.

This script is **not** part of the installed package and is not needed at
run time. Re-run it only to refresh the vendored catalogue. It needs network
access and ``requests`` (installed with ScopeSim, but not a direct dependency
of ScopeSim_Templates).

Ported from ``fetch_catalogue`` in the MAVIS 2026 workshop ``hudf.py``.

Usage (from the repository root)::

    python scripts/make_hudf_catalogue.py               # 2.5 arcmin cone
    python scripts/make_hudf_catalogue.py --radius 1.5 --out my_cat.fits
"""

import argparse
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table


#: HUDF ACS field centre (Beckwith et al. 2006, AJ 132, 1729), degrees.
HUDF_CENTRE = (53.1625, -27.79142)
VIZIER_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
OUT_PATH = (Path(__file__).resolve().parents[1] / "scopesim_templates"
            / "extragalactic" / "data" / "hudf_catalogue.fits")


def _asu_query(source, columns, centre, radius_arcmin, constraints=None,
               max_rows=200000, timeout=180):
    """Fetch one VizieR table as a TSV cone search, and parse it."""
    import requests

    params = {
        "-source": source,
        "-c": f"{centre[0]} {centre[1]}",
        "-c.rm": f"{radius_arcmin}",
        "-out.max": str(max_rows),
        "-out": ",".join(columns),
    }
    params.update(constraints or {})
    reply = requests.get(VIZIER_ASU, params=params, timeout=timeout)
    reply.raise_for_status()
    if "Error=" in reply.text:
        bad = next(ln for ln in reply.text.splitlines() if "Error=" in ln)
        raise RuntimeError(f"VizieR refused {source}: {bad.strip()}")

    rows = []
    for line in reply.text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < len(columns):
            continue
        try:                       # skip the unit and dashed separator rows
            float(parts[0])
        except ValueError:
            continue
        rows.append(parts)

    if not rows:
        raise RuntimeError(f"VizieR returned no rows for {source}")

    table = Table()
    for i, name in enumerate(columns):
        values = [row[i] for row in rows]
        try:
            table[name] = np.array(
                [float(v) if v not in ("", "-") else np.nan for v in values])
        except ValueError:
            table[name] = np.array(values)
    return table


def build_catalogue(radius_arcmin=2.5, centre=HUDF_CENTRE, match_radius=0.5,
                    out_path=OUT_PATH):
    """Query, cross-match and write the catalogue. Returns the Table."""
    print(f"Querying van der Wel+2012 within {radius_arcmin}' ...")
    structure = _asu_query(
        "J/ApJS/203/24/table2",
        ["RAJ2000", "DEJ2000", "ID", "Q", "mag", "r", "n", "q", "PA"],
        centre, radius_arcmin, constraints={"F": "H", "Fld": "G"})
    print(f"  {len(structure)} rows")
    print("Querying Skelton+2014 ...")
    photometry = _asu_query(
        "J/ApJS/214/24/3dhstall",
        ["RAJ2000", "DEJ2000", "zsp", "zpk", "F606W", "F160W", "logM", "Use"],
        centre, radius_arcmin, constraints={"Field": "GOODS-S"})
    print(f"  {len(photometry)} rows")

    good = ((structure["Q"] == 0) & np.isfinite(structure["mag"])
            & np.isfinite(structure["r"]) & np.isfinite(structure["n"])
            & np.isfinite(structure["q"]) & (structure["r"] > 0))
    structure = structure[good]

    left = SkyCoord(structure["RAJ2000"] * u.deg, structure["DEJ2000"] * u.deg)
    right = SkyCoord(photometry["RAJ2000"] * u.deg,
                     photometry["DEJ2000"] * u.deg)
    idx, sep2d, _ = left.match_to_catalog_sky(right)
    matched = sep2d.arcsec < match_radius
    structure, partner = structure[matched], photometry[idx[matched]]
    sep = sep2d.arcsec[matched]

    z_spec = np.isfinite(partner["zsp"]) & (partner["zsp"] > 0)
    redshift = np.where(z_spec, partner["zsp"], partner["zpk"])

    with np.errstate(divide="ignore", invalid="ignore"):
        # 3D-HST fluxes are on an AB zero point of 25
        m_f606w = 25.0 - 2.5 * np.log10(
            np.where(partner["F606W"] > 0, partner["F606W"], np.nan))

    out = Table({
        "id": np.asarray(structure["ID"]).astype(np.int32),
        "ra": np.asarray(structure["RAJ2000"]),
        "dec": np.asarray(structure["DEJ2000"]),
        "z": redshift,
        "z_spec": z_spec,
        "m_F160W": np.asarray(structure["mag"]),
        "r_eff_arcsec": np.asarray(structure["r"]),
        "sersic_n": np.asarray(structure["n"]),
        "axis_ratio": np.asarray(structure["q"]),
        "position_angle_deg": np.asarray(structure["PA"]),
        "log_mstar": np.asarray(partner["logM"]),
        "m_F606W": m_f606w,
        "match_sep_arcsec": sep,
    })
    out = out[np.isfinite(out["z"]) & (out["z"] > 0)]
    out.sort("m_F160W")

    for col, unit in [("ra", u.deg), ("dec", u.deg), ("m_F160W", u.mag),
                      ("r_eff_arcsec", u.arcsec), ("m_F606W", u.mag),
                      ("position_angle_deg", u.deg),
                      ("match_sep_arcsec", u.arcsec)]:
        out[col].unit = unit

    out.meta.update({
        "CENTRA": centre[0],         # deg, Beckwith+2006
        "CENTDEC": centre[1],        # deg
        "RADIUS": radius_arcmin,     # arcmin
        "MATCHRAD": match_radius,    # arcsec
        "STRUCT": "van der Wel+2012 ApJS 203 24, J/ApJS/203/24 (F160W, Q=0)",
        "PHOTOM": "Skelton+2014 ApJS 214 24, J/ApJS/214/24 (GOODS-S)",
        "PACONV": "position_angle_deg: 0 = North, 90 = East",
        "MAGSYS": "AB",
    })

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.write(out_path, overwrite=True)
    print(f"Wrote {len(out)} galaxies to {out_path} "
          f"({out_path.stat().st_size / 1024:.1f} kB)")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--radius", type=float, default=2.5,
                        help="cone radius in arcmin (default 2.5)")
    parser.add_argument("--match-radius", type=float, default=0.5,
                        help="cross-match radius in arcsec (default 0.5)")
    parser.add_argument("--out", type=Path, default=OUT_PATH,
                        help="output FITS path (default: vendored location)")
    args = parser.parse_args()
    build_catalogue(radius_arcmin=args.radius, match_radius=args.match_radius,
                    out_path=args.out)


if __name__ == "__main__":
    main()
