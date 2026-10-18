# -*- coding: utf-8 -*-
"""A hypothetical unseen stellar population on random orbits about Sgr A*.

Fills in what the S-star catalogue
(:func:`~scopesim_templates.stellar.galactic_centre`) cannot: the *unseen*
stars. A Kroupa (2001) IMF is sampled, every star is given a random
Keplerian orbit about the black hole, and the orbits are evaluated at the
requested epoch -- so a frame shows plausible background crowding and
plausible, position-correlated proper motions instead of an empty sky.

How much mass to put there
--------------------------
There is **no positive measurement** of the stellar mass inside the central
arcsecond, only an upper limit. GRAVITY Collaboration 2022 (A&A 657, L12;
arXiv:2112.07478) constrain any extended mass inside S2's apocentre,
r = 0.23 arcsec, to **<~3000 Msun (1 sigma)** for plausible density
profiles, and only detect a tentative excess at r >~ 2.5 arcsec. Anything
from ~0 to a few thousand Msun is therefore allowed. The default here,
100 Msun inside r = 1 arcsec, is deliberately conservative; pass
``total_mass=GRAVITY_2022_MASS_LIMIT`` (or ``"gravity2022"``) for the
3000 Msun ceiling. Note that the template spreads ``total_mass`` over
``r < r_max`` (default 1 arcsec), whereas the GRAVITY limit applies inside
0.23 arcsec, so 3000 Msun inside 1 arcsec is a generous bracket rather than
the literal limit.

What the population looks like
------------------------------
Measured in the workshop with a 384 mas field and an m_K < 24 cut:

===============  =========  =========  ==============  ===================
total mass       sampled    remnants   in the field    proper motion
===============  =========  =========  ==============  ===================
100 Msun            180        14           10         12-67 mas/yr
3000 Msun          5391       413          312          1-137 mas/yr
===============  =========  =========  ==============  ===================

In both cases the brightest background star is m_K ~ 19.7 (fainter than
S301) with a median near m_K 22.6: the turn-off caps the top of the main
sequence, so mass is spread over many faint stars.

Orbits
------
Elements are drawn as an isotropic, relaxed cluster around a point mass:

* semi-major axis from a cusp ``rho ~ r^-gamma``, i.e. ``N(<a) ~
  a^(3-gamma)``, gamma = 1.4 (Schoedel et al. 2018, Gallego-Cano et al.
  2018);
* thermal eccentricity, ``f(e) = 2e``; isotropic orientation;
* uniform mean anomaly at :data:`REFERENCE_EPOCH`; period from Kepler's
  third law with M_BH = 4.297e6 Msun (GRAVITY 2022).

The orbital elements are drawn once per ``seed``, so calling
:func:`galactic_centre_background` with the same ``seed`` and different
``time`` values gives successive epochs of the *same* population.

Assumptions
-----------
* Semi-major axes are sampled from the *space density* profile; for
  eccentric orbits this puts slightly too many stars at small radii.
* Only main-sequence stars are rendered. Masses above ``turnoff_mass``
  (1.2 Msun, roughly a 5 Gyr population) are dropped as dark remnants.
* Magnitudes come from the zero-age, solar-metallicity Mamajek dwarf
  sequence (2MASS Ks), with a single A_K = 2.42 (Fritz et al. 2011) and
  R0 = 8277 pc.
* Spectra are six Pickles dwarf templates binned by mass and normalised to
  0 mag in ``filter_name``; each star's flux is carried by its ``weight``.
  Unlike :func:`galactic_centre`, spectra are *not* RV-shifted per star (that
  would need one spectrum per star); the line-of-sight velocity is stored in
  the ``rv`` column instead.
* No mass segregation, binaries, remnants, or extinction variation.

Ported from ``gc_population.py`` of the MAVIS 2026 workshop feasibility
studies, with the ``scopesim_targets`` IMF and stellar-parameter lookup
replaced by self-contained code (:mod:`._population`).
"""

import numpy as np
from astropy import units as u
from astropy.table import Table

import pyckles
from spextra import Spextrum

from ..rc import Source
from ..utils.general_utils import add_function_call_str
from ._orbits import (sky_state, to_jyear, M_BH_SGR_A, R0_SGR_A,
                      RA_SGR_A_STAR, DEC_SGR_A_STAR)
from ._population import (sample_broken_powerlaw, broken_powerlaw_mean,
                          dwarf_absmag, seed_for_meta)


__all__ = ["galactic_centre_background", "GRAVITY_2022_MASS_LIMIT"]


#: [Msun] 1-sigma upper limit on extended mass inside S2's apocentre
#: (GRAVITY Collaboration 2022, A&A 657, L12)
GRAVITY_2022_MASS_LIMIT = 3000.0

#: K-band extinction toward the Galactic Centre (Fritz et al. 2011)
A_K = 2.42
#: [Msun] main-sequence turn-off of an old (~5 Gyr) population
TURNOFF_MASS = 1.2
#: slope of the stellar cusp, rho ~ r^-gamma
GAMMA = 1.4
#: [Msun] IMF limits: hydrogen-burning limit to the IMF top
MIN_MASS, IMF_TOP = 0.08, 60.0
#: [Julian yr] epoch at which mean anomalies are uniform; fixed so that a
#: given seed describes the same population at every epoch
REFERENCE_EPOCH = 2030.0

#: Pickles dwarf templates spanning the main sequence below the turn-off,
#: as (lower mass bound [Msun], template name), in increasing mass.
TEMPLATE_BINS = [(0.00, "M5V"), (0.12, "M3V"), (0.30, "M0V"),
                 (0.50, "K5V"), (0.70, "G5V"), (0.95, "F8V")]


def _parse_total_mass(total_mass):
    if isinstance(total_mass, str):
        if total_mass.strip().lower() in ("gravity2022", "upper_limit"):
            return GRAVITY_2022_MASS_LIMIT
        raise ValueError(f"Unknown total_mass alias {total_mass!r}; use a "
                         "number [Msun] or 'gravity2022'.")
    return float(total_mass)


def _sample_population(total_mass, r_max, gamma, turnoff_mass, rng):
    """Draw masses and orbital elements. Same ``rng`` state -> same stars."""
    mean_mass = broken_powerlaw_mean(MIN_MASS, IMF_TOP)
    n_sampled = int(round(total_mass / mean_mass))
    mass = sample_broken_powerlaw(n_sampled, MIN_MASS, IMF_TOP, rng)

    # Draw every element for every sampled star *before* dropping remnants,
    # so the stream of random numbers does not depend on turnoff_mass.
    radius = r_max * rng.random(n_sampled) ** (1.0 / (3.0 - gamma))  # arcsec
    ecc = np.sqrt(rng.random(n_sampled))                          # f(e) = 2e
    incl = np.degrees(np.arccos(rng.uniform(-1, 1, n_sampled)))
    node = rng.uniform(0, 360, n_sampled)
    argp = rng.uniform(0, 360, n_sampled)
    phase = rng.random(n_sampled)

    a_au = radius * R0_SGR_A               # arcsec * pc = au
    period = np.sqrt(a_au**3 / M_BH_SGR_A)  # yr
    alive = mass <= turnoff_mass
    return {
        "mass": mass[alive], "a_mas": radius[alive] * 1e3, "e": ecc[alive],
        "incl": incl[alive], "node": node[alive], "argp": argp[alive],
        "t_peri": (REFERENCE_EPOCH - phase * period)[alive],
        "period": period[alive],
        "n_sampled": n_sampled, "n_remnant": int((~alive).sum()),
        "mean_mass": mean_mass,
    }


def _background_catalogue(time, total_mass=100.0, r_max=1.0, gamma=GAMMA,
                          turnoff_mass=TURNOFF_MASS, extinction_k=A_K,
                          mag_limit=None, seed=None):
    """Positions, velocities and magnitudes at ``time`` (no spectra).

    Returns an ``astropy.table.Table`` with columns ``x``, ``y`` [arcsec],
    ``m_k`` [mag], ``mass`` [Msun], ``rv``, ``pm_x``, ``pm_y`` and the
    orbital elements; bookkeeping numbers are in ``table.meta``.
    """
    total_mass = _parse_total_mass(total_mass)
    if total_mass <= 0:
        raise ValueError("total_mass must be positive.")
    rng = np.random.default_rng(seed)
    pop = _sample_population(total_mass, float(r_max), float(gamma),
                             float(turnoff_mass), rng)
    jyear = to_jyear(time)

    distmod = 5 * np.log10(R0_SGR_A / 10.0) + extinction_k
    m_k = dwarf_absmag(pop["mass"], "M_Ks") + distmod

    elements = [pop[k] for k in ("a_mas", "e", "incl", "node", "argp",
                                 "t_peri", "period")]
    now = sky_state(*elements, jyear)
    # proper motion from a symmetric finite difference over +/- 1 day
    step = 1.0 / 365.25
    ahead = sky_state(*elements, jyear + step)
    behind = sky_state(*elements, jyear - step)
    pm_x = (ahead["east"] - behind["east"]) / (2 * step)
    pm_y = (ahead["north"] - behind["north"]) / (2 * step)

    tbl = Table({
        "x": now["east"] * 1e-3 * u.arcsec,
        "y": now["north"] * 1e-3 * u.arcsec,
        "m_k": m_k * u.mag,
        "mass": pop["mass"] * u.Msun,
        "rv": now["v_r"] * u.km / u.s,
        "pm_x": pm_x * u.mas / u.yr,
        "pm_y": pm_y * u.mas / u.yr,
        "a": pop["a_mas"] * u.mas,
        "e": pop["e"],
        "incl": pop["incl"] * u.deg,
        "node": pop["node"] * u.deg,
        "argp": pop["argp"] * u.deg,
        "t_peri": pop["t_peri"],
        "period": pop["period"] * u.yr,
    })
    tbl["star_id"] = np.arange(len(tbl))
    if mag_limit is not None:
        tbl = tbl[tbl["m_k"].value <= mag_limit]

    tbl.meta.update({
        "total_mass": total_mass,
        "r_max": float(r_max),
        "gamma": float(gamma),
        "turnoff_mass": float(turnoff_mass),
        "extinction_k": float(extinction_k),
        "n_sampled": pop["n_sampled"],
        "n_remnant": pop["n_remnant"],
        "mean_mass": pop["mean_mass"],
        "epoch_jyear": jyear,
        "reference_epoch": REFERENCE_EPOCH,
    })
    return tbl


def _template_index(mass):
    """Index into :data:`TEMPLATE_BINS` for each mass."""
    edges = np.array([lo for lo, _ in TEMPLATE_BINS])
    return np.searchsorted(edges, np.asarray(mass), side="right") - 1


@add_function_call_str
def galactic_centre_background(time,
                               total_mass=100.0,
                               r_max=1.0,
                               gamma=GAMMA,
                               turnoff_mass=TURNOFF_MASS,
                               extinction_k=A_K,
                               mag_limit=None,
                               filter_name="Paranal/NACO.Ks",
                               seed=None):
    """
    Source of a hypothetical unseen main-sequence population around Sgr A*.

    Kroupa-IMF stars on random Keplerian orbits about Sgr A*, evaluated at
    ``time``. See the module docstring for the physical assumptions and
    the GRAVITY 2022 mass constraint.

    Parameters
    ----------
    time : astropy.time.Time, datetime, str, or float MJD
        Epoch. Same conventions as
        :func:`~scopesim_templates.stellar.galactic_centre`.
    total_mass : float or str
        [Msun] Total *initial* stellar mass inside ``r_max`` (remnants
        above the turn-off are counted in it but not rendered). Default 100.
        ``"gravity2022"`` (or :data:`GRAVITY_2022_MASS_LIMIT`) gives the
        3000 Msun upper limit of GRAVITY Collaboration 2022.
    r_max : float
        [arcsec] Outer semi-major axis of the population.
    gamma : float
        Slope of the density cusp, rho ~ r^-gamma.
    turnoff_mass : float
        [Msun] Stars above this are dropped as dark remnants.
    extinction_k : float
        [mag] K-band extinction applied to every star.
    mag_limit : float, optional
        Drop stars fainter than this apparent K magnitude. Applied after
        sampling, so it does not change which stars are drawn.
    filter_name : str
        Filter in which the magnitudes are applied. Default
        ``"Paranal/NACO.Ks"``, as in :func:`galactic_centre`.
    seed : int, np.random.Generator, or None
        Random seed. ``None`` (default) gives a different population on
        every call. Use the same integer seed at different ``time`` values
        to follow one population through time.

    Returns
    -------
    scopesim.Source
        Columns ``x``, ``y`` [arcsec offsets from Sgr A*, x = East],
        ``ref``, ``weight``, ``m_k``, ``mass``, ``rv`` [km/s], ``pm_x``,
        ``pm_y`` [mas/yr], ``star_id`` and the orbital elements.

    See Also
    --------
    galactic_centre : the S-star field, which can include this population
        via ``background_mass``.

    Examples
    --------
    >>> from scopesim_templates.stellar import galactic_centre_background
    >>> src = galactic_centre_background("2031-10-22", total_mass=3000,
    ...                                  mag_limit=24, seed=42)
    """
    tbl = _background_catalogue(time, total_mass=total_mass, r_max=r_max,
                                gamma=gamma, turnoff_mass=turnoff_mass,
                                extinction_k=extinction_k,
                                mag_limit=mag_limit, seed=seed)
    if not len(tbl):
        raise ValueError("No background stars were generated; increase "
                         "total_mass or relax mag_limit.")

    # One spectrum per occupied template bin, normalised to 0 mag
    idx = _template_index(tbl["mass"].value)
    used = np.unique(idx)
    pickles_lib = pyckles.SpectralLibrary("pickles", return_style="synphot")
    spectra = [Spextrum(modelclass=pickles_lib[TEMPLATE_BINS[i][1]])
               .scale_to_magnitude(filter_curve=filter_name,
                                   amplitude=0 * u.mag)
               for i in used]
    tbl["ref"] = np.searchsorted(used, idx)
    tbl["weight"] = 10 ** (-0.4 * tbl["m_k"].value)
    tbl["spec_types"] = [TEMPLATE_BINS[i][1] for i in idx]

    meta = dict(tbl.meta)
    tbl.meta.clear()
    src = Source(spectra=spectra, table=tbl)
    src.meta.update(meta)
    src.meta.update({
        "description": "Hypothetical background population around Sgr A*",
        "epoch": str(time),
        "seed": seed_for_meta(seed),
        "ra": RA_SGR_A_STAR,
        "dec": DEC_SGR_A_STAR,
        "m_bh": M_BH_SGR_A,
        "r0_pc": R0_SGR_A,
        "filter_name": filter_name,
    })
    return src
