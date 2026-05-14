"""Evaluate Gillessen+2017 Galactic Centre stellar orbits at arbitrary epochs.

Loads the catalogue (CDS J/ApJ/837/30 table 3, 40 stars orbiting Sgr A*) and
returns, for any time, each star's projected sky position (d_ra, d_dec in
arcsec, east-positive / north-positive) and line-of-sight radial velocity
(km/s, positive = receding), alongside catalogue Kmag and spectral type.

Public entry point: ``stars_at_time(time)``.

Vendored from the standalone ``gillessen`` script collection. Catalogue
reference: Gillessen et al. 2017, ApJ 837, 30 — CDS table J/ApJ/837/30.
The companion FITS catalogue ``J_ApJ_837_30_table3.dat.fits`` is shipped
alongside this module.
"""

from __future__ import annotations

import warnings
from datetime import datetime
from numbers import Real
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.table import QTable, Table
from astropy.time import Time

DEFAULT_CATALOGUE = Path(__file__).resolve().parent / "J_ApJ_837_30_table3.dat.fits"
DEFAULT_R0 = 8.32 * u.kpc  # Gillessen+2017


def _load_catalogue(path: Path | str = DEFAULT_CATALOGUE) -> Table:
    return Table.read(path)


def _to_jyear(t) -> float:
    """Convert input time to Julian year (J2000-based, matches catalogue Tp).

    Accepts ``astropy.time.Time``, ``datetime``, ISO string, or numeric MJD.
    """
    if isinstance(t, Time):
        return float(t.jyear)
    if isinstance(t, datetime):
        return float(Time(t).jyear)
    if isinstance(t, Real) and not isinstance(t, bool):
        return float(Time(float(t), format="mjd").jyear)
    return float(Time(t).jyear)


def _solve_kepler(M, e, tol: float = 1e-12, max_iter: int = 50) -> np.ndarray:
    """Solve Kepler's equation ``M = E - e sin E`` for ``E``.

    Vectorised Newton iteration; ``M`` is wrapped to ``[-pi, pi)``.
    """
    M = np.asarray(M, dtype=float)
    e = np.asarray(e, dtype=float)
    M = np.mod(M + np.pi, 2.0 * np.pi) - np.pi
    E = M + e * np.sin(M)
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = f / fp
        E = E - dE
        if np.all(np.abs(dE) < tol):
            break
    return E


def _orbital_state(a, e, n, E):
    """Position and velocity in the orbital plane.

    Pericenter on +x. Returns ``(x_orb, y_orb, vx_orb, vy_orb)`` in
    ``(a-units, a-units / time-units-of-n)``.
    """
    cosE = np.cos(E)
    sinE = np.sin(E)
    sqrt1me2 = np.sqrt(1.0 - e * e)
    x_orb = a * (cosE - e)
    y_orb = a * sqrt1me2 * sinE
    Edot = n / (1.0 - e * cosE)
    vx_orb = -a * sinE * Edot
    vy_orb = a * sqrt1me2 * cosE * Edot
    return x_orb, y_orb, vx_orb, vy_orb


def _thiele_innes(Omega, inc, omega):
    """Thiele-Innes constants (A, B, F, G, C, H), all angles in radians.

    Mapping (Murray & Dermott 1999):
        d_dec (north) = A*x + F*y
        d_ra  (east)  = B*x + G*y
        z (away from observer) = C*x + H*y
    """
    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(omega), np.sin(omega)
    A = cO * cw - sO * sw * ci
    B = sO * cw + cO * sw * ci
    F = -cO * sw - sO * cw * ci
    G = -sO * sw + cO * cw * ci
    C = sw * si
    H = cw * si
    return A, B, F, G, C, H


def _arcsec_per_yr_to_km_per_s(v_arcsec_per_yr: np.ndarray, R0: u.Quantity) -> u.Quantity:
    """Convert angular speed at distance ``R0`` to linear km/s."""
    return (v_arcsec_per_yr * (u.arcsec / u.yr) * R0).to(
        u.km / u.s, equivalencies=u.dimensionless_angles()
    )


def stars_at_time(time, catalogue: Path | str | None = None, R0: u.Quantity = DEFAULT_R0) -> QTable:
    """Predict positions and radial velocities of all catalogue stars at ``time``.

    Parameters
    ----------
    time : astropy.time.Time, datetime, ISO string, or float MJD
        Epoch at which to evaluate the orbits.
    catalogue : path-like, optional
        Path to ``J_ApJ_837_30_table3.dat.fits``. Defaults to the file
        shipped next to this module.
    R0 : astropy.units.Quantity, optional
        Distance to Sgr A*. Default 8.32 kpc (Gillessen+2017). Used only to
        convert orbital angular velocity to line-of-sight km/s.

    Returns
    -------
    astropy.table.QTable
        Columns: ``Star`` (str), ``d_ra`` [arcsec, east-positive],
        ``d_dec`` [arcsec, north-positive], ``RV`` [km/s, positive=receding],
        ``Kmag`` [mag], ``SpT`` (str: ``e`` early / ``l`` late).
        Hyperbolic / open-orbit rows (e.g. S111) are skipped with a warning.
    """
    cat = _load_catalogue(catalogue if catalogue is not None else DEFAULT_CATALOGUE)
    t_jyear = _to_jyear(time)

    a = np.asarray(cat["a"], dtype=float)
    e = np.asarray(cat["e"], dtype=float)
    per = np.asarray(cat["Per"], dtype=float)
    valid = (a > 0) & (e >= 0) & (e < 1) & (per > 0)
    if not valid.all():
        skipped = [
            (s.decode().strip() if isinstance(s, bytes) else str(s).strip())
            for s in np.asarray(cat["Star"])[~valid]
        ]
        warnings.warn(
            f"Skipping hyperbolic / open-orbit rows: {skipped}",
            UserWarning,
            stacklevel=2,
        )

    cat = cat[valid]
    a = a[valid]
    e = e[valid]
    per = per[valid]
    inc = np.deg2rad(np.asarray(cat["i"], dtype=float))
    Omega = np.deg2rad(np.asarray(cat["Omega"], dtype=float))
    omega = np.deg2rad(np.asarray(cat["w"], dtype=float))
    Tp = np.asarray(cat["Tp"], dtype=float)

    n = 2.0 * np.pi / per  # rad / Julian yr
    M = n * (t_jyear - Tp)
    E = _solve_kepler(M, e)
    x_orb, y_orb, vx_orb, vy_orb = _orbital_state(a, e, n, E)
    A, B, F, G, C, H = _thiele_innes(Omega, inc, omega)

    d_dec = A * x_orb + F * y_orb
    d_ra = B * x_orb + G * y_orb
    vz_arcsec_per_yr = C * vx_orb + H * vy_orb
    rv = _arcsec_per_yr_to_km_per_s(vz_arcsec_per_yr, R0)

    def _decode(x):
        return x.decode().strip() if isinstance(x, bytes) else str(x).strip()

    star_names = [_decode(s) for s in np.asarray(cat["Star"])]
    spt = [_decode(s) for s in np.asarray(cat["SpT"])]

    out = QTable()
    out["Star"] = star_names
    out["d_ra"] = d_ra * u.arcsec
    out["d_dec"] = d_dec * u.arcsec
    out["RV"] = rv
    out["Kmag"] = np.asarray(cat["Kmag"], dtype=float) * u.mag
    out["SpT"] = spt
    return out
