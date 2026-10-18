# -*- coding: utf-8 -*-
"""Generator-based stellar-population helpers for the science-case templates.

Self-contained replacements for the pieces of ``scopesim_targets`` used by
the MAVIS 2026 workshop code (``BrokenPowerlaw`` IMF and
``StellarParameters`` lookup), so that no new dependency is needed:

* :func:`sample_broken_powerlaw` -- exact inverse-CDF sampling of a
  continuous broken power-law IMF from an ``np.random.Generator``. Unlike
  :mod:`scopesim_templates.stellar.imf` it never touches the global NumPy
  RNG.
* :func:`broken_powerlaw_mean` -- the (correctly normalised) mean mass.
* :func:`dwarf_absmag` -- absolute magnitude of a zero-age, solar-metallicity
  dwarf of given mass, interpolated in the Mamajek table already shipped
  with the package (``mamajek_alt.dat``; Pecaut & Mamajek 2013, ApJS 208, 9,
  version 2019.3.22).
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
from astropy.table import Table


__all__ = ["KROUPA_2001", "sample_broken_powerlaw", "broken_powerlaw_mean",
           "dwarf_absmag", "dwarf_table", "seed_for_meta"]


DIRNAME = Path(__file__).parent

#: Kroupa (2001) IMF, dN/dm ~ m^-alpha: (mass breaks [Msun], slopes alpha).
#: Identical to the parameters used by the workshop code.
KROUPA_2001 = ((0.01, 0.08, 0.5, 150.0), (0.3, 1.3, 2.3))


def seed_for_meta(seed):
    """A FITS/YAML-friendly representation of ``seed`` for ``Source.meta``."""
    if seed is None or isinstance(seed, (int, np.integer)):
        return None if seed is None else int(seed)
    return repr(seed)


def _segments(m_min, m_max, imf):
    """Clip a broken power law to [m_min, m_max].

    Returns segment edges (lo, hi), slopes, and the relative (unnormalised,
    but continuous across the breaks) amplitude of each segment.
    """
    breaks, slopes = (np.asarray(v, dtype=float) for v in imf)
    if not (breaks[0] <= m_min < m_max <= breaks[-1]):
        raise ValueError(f"Mass range ({m_min}, {m_max}) must lie inside the "
                         f"IMF range ({breaks[0]}, {breaks[-1]}).")
    # amplitudes that make xi(m) continuous: A_k m_k^-a_k = A_{k-1} m_k^-a_{k-1}
    amps = np.ones(len(slopes))
    for k in range(1, len(slopes)):
        amps[k] = amps[k - 1] * breaks[k] ** (slopes[k] - slopes[k - 1])
    lo = np.clip(breaks[:-1], m_min, m_max)
    hi = np.clip(breaks[1:], m_min, m_max)
    keep = hi > lo
    return lo[keep], hi[keep], slopes[keep], amps[keep]


def _power_integral(lo, hi, power):
    """Integral of m^power from lo to hi, element-wise."""
    out = np.empty_like(lo)
    log = np.isclose(power, -1.0)
    out[log] = np.log(hi[log] / lo[log])
    p1 = power[~log] + 1
    out[~log] = (hi[~log] ** p1 - lo[~log] ** p1) / p1
    return out


def broken_powerlaw_mean(m_min, m_max, imf=KROUPA_2001):
    """Mean stellar mass [Msun] of a broken power-law IMF on [m_min, m_max].

    Computed as int(m xi) / int(xi), i.e. correctly normalised also for a
    truncated IMF (see the ``scopesim_targets`` gotcha about ``expect()``).
    """
    lo, hi, slopes, amps = _segments(m_min, m_max, imf)
    number = np.sum(amps * _power_integral(lo, hi, -slopes))
    mass = np.sum(amps * _power_integral(lo, hi, 1.0 - slopes))
    return float(mass / number)


def sample_broken_powerlaw(n, m_min, m_max, rng, imf=KROUPA_2001):
    """Draw ``n`` masses [Msun] from a broken power-law IMF.

    Exact inverse-CDF sampling: a segment is chosen with probability equal
    to its share of the number integral, then the mass is drawn from the
    analytic inverse CDF of ``m^-alpha`` within that segment.

    Parameters
    ----------
    n : int
    m_min, m_max : float
        Mass range [Msun]; must lie inside the IMF's own range.
    rng : np.random.Generator
    imf : tuple
        ``(breaks, slopes)`` with ``len(breaks) == len(slopes) + 1``.
        Default :data:`KROUPA_2001`.

    Returns
    -------
    ndarray of float
    """
    n = int(n)
    if n <= 0:
        return np.empty(0)
    lo, hi, slopes, amps = _segments(m_min, m_max, imf)
    weights = amps * _power_integral(lo, hi, -slopes)
    seg = rng.choice(len(lo), size=n, p=weights / weights.sum())
    u = rng.random(n)

    lo, hi, a = lo[seg], hi[seg], slopes[seg]
    out = np.empty(n)
    log = np.isclose(a, 1.0)
    out[log] = lo[log] * (hi[log] / lo[log]) ** u[log]
    p1 = 1.0 - a[~log]
    lo_p, hi_p = lo[~log] ** p1, hi[~log] ** p1
    out[~log] = (lo_p + u[~log] * (hi_p - lo_p)) ** (1.0 / p1)
    return out


@lru_cache(maxsize=1)
def dwarf_table():
    """Mamajek dwarf sequence (SpT, Mv, M_J, M_Ks, Msun), sorted by mass."""
    tbl = Table.read(DIRNAME / "mamajek_alt.dat", format="ascii.fixed_width")
    tbl.sort("Msun")
    return tbl


def dwarf_absmag(mass, band="Mv"):
    """Absolute magnitude of a ZAMS dwarf of ``mass`` [Msun].

    Linear interpolation in the Mamajek table. ``band`` is one of ``"Mv"``
    (Johnson V), ``"M_J"`` or ``"M_Ks"`` (2MASS). Outside the tabulated
    range (0.01 - ~20 Msun for the IR bands) the end values are used.
    """
    tbl = dwarf_table()
    col = tbl[band]
    good = ~np.ma.getmaskarray(col)
    return np.interp(np.asarray(mass, dtype=float),
                     np.asarray(tbl["Msun"], dtype=float)[good],
                     np.asarray(col, dtype=float)[good])
