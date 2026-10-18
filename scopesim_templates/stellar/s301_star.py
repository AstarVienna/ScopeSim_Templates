# -*- coding: utf-8 -*-
"""S301 -- the star that feels the spin of Sgr A* -- and the Sgr A* NIR source.

Orbit model and ``Source`` builder for S301, the 8.7-year star discovered
with GRAVITY+ (Mang et al. 2026, "Discovery of a star sensitive to the spin
of Sgr A*", Nature; arXiv:2607.12664).

Two orbit solutions
-------------------
Astrometry alone cannot tell a stellar orbit from its mirror image through
the plane of the sky: rotating both the node and the argument of periapsis
by 180 deg leaves the *projected* ellipse untouched but flips the sign of the
line-of-sight velocity. The discovery paper therefore quotes two equally
good solutions (Extended Data Table 2), available here as
``S301_SOLUTIONS[1]`` and ``S301_SOLUTIONS[2]``. Near the 2031.81
pericentre they differ by ~17500 km/s in radial velocity, i.e. Brackett
gamma lands at either ~2.24 um or ~2.12 um: one spectrum breaks the
degeneracy.

The pericentre crossing takes only ~19 hours, and the radial velocity
sweeps at up to 13500 km/s per day there. For spectroscopy the useful epoch
is the *stationary point* of the radial-velocity curve, about one day after
pericentre (see the workshop documentation).

Validation (from the elements, with M_BH and R0 of GRAVITY 2022):
a = 687 au; Kepler period 8.69 yr (fitted 8.68 +/- 0.11); pericentre
136 / 143 R_s and peak speed 25600 / 25000 km/s for solutions 1 / 2, as in
the paper.

What is *not* modelled
----------------------
* No Lense-Thirring term (the actual science goal), no Schwarzschild
  precession in the propagation (2.0 deg/orbit,
  ``KeplerOrbit.precession_per_orbit``), no light bending, Shapiro or
  Roemer delay.
* The redshift is the leading-order sum of Doppler, transverse-Doppler and
  gravitational terms (``relativistic=True``), not an integrated null
  geodesic. At pericentre the relativistic terms add ~2200 km/s, so they
  are on by default.
* Sgr A*'s observed m_K = 17.2 is a *derived* median (GRAVITY median
  dereddened 2.2 um flux ~0.8 mJy, reddened by A_K = 2.42); it flares to
  m_K ~ 15 several times a night.
* The stellar template is the BOSZ F0V high-resolution model (the paper's
  best guess is F1.5V, which BOSZ does not have). Its Brackett gamma line
  (EW 15.8 A, 32 per cent core depth) is the deepest feature in the K band.

References
----------
Mang et al. 2026, Nature (arXiv:2607.12664) -- orbital elements (Extended
    Data Table 2), m_K = 19.3 +/- 0.3, "late A-type or early F-type".
GRAVITY Collaboration 2022, A&A 657, L12 -- M_BH = 4.297e6 Msun,
    R0 = 8277 pc.

Ported from ``s301.py`` of the MAVIS 2026 workshop feasibility studies.
"""

import numpy as np
from astropy import units as u
from astropy.table import Table

from spextra import Spextrum

from ..rc import Source
from ..utils.general_utils import add_function_call_str
from ._orbits import (KeplerOrbit, to_jyear, R0_SGR_A, M_BH_SGR_A,
                      RA_SGR_A_STAR, DEC_SGR_A_STAR)


__all__ = ["s301", "s301_orbit", "s301_spectrum", "sgr_a_star_spectrum",
           "S301_SOLUTIONS"]


#: Solution 1 of Mang et al. 2026, Extended Data Table 2.
S301_SOLUTION_1 = KeplerOrbit(label="S301 solution 1", a=83.0, e=0.9832,
                              incl=124.09, node=73.8, argp=293.4,
                              t_peri=2023.126, period=8.68)

#: Solution 2 -- the (near) mirror image, equally consistent with astrometry.
S301_SOLUTION_2 = KeplerOrbit(label="S301 solution 2", a=83.0, e=0.9824,
                              incl=122.84, node=256.9, argp=115.1,
                              t_peri=2023.125, period=8.68)

#: Both solutions, keyed by number.
S301_SOLUTIONS = {1: S301_SOLUTION_1, 2: S301_SOLUTION_2}

#: 1-sigma uncertainties of solution 1 (same table), for error budgets.
S301_SOLUTION_1_ERRORS = dict(a=0.7, e=0.0010, incl=1.10, node=3.5,
                              argp=2.2, t_peri=0.010, period=0.11)

#: m_K = 19.3 +/- 0.3 (Mang et al. 2026)
S301_MK = 19.3
#: Closest BOSZ high-resolution template to the paper's F1.5V
S301_TEMPLATE = "bosz/hr/f0v"

#: Median observed K magnitude of Sgr A*'s NIR counterpart (derived, see
#: module docstring) and its featureless power-law slope, F_nu ~ nu^alpha.
SGR_A_STAR_MK = 17.2
SGR_A_STAR_ALPHA = -0.6

#: Brackett gamma (7-4) and delta (8-4), vacuum wavelengths
BRACKETT_GAMMA = 2.16612 * u.um
BRACKETT_DELTA = 1.94509 * u.um


def s301_orbit(solution=1):
    """Return one of the two published S301 orbit solutions.

    Parameters
    ----------
    solution : {1, 2}

    Returns
    -------
    KeplerOrbit
        Frozen dataclass with ``sky_position(jyear)`` (east, north in mas),
        ``radial_velocity(jyear)``, ``apparent_radial_velocity(jyear)``
        [km/s], ``next_pericentre(after)``, ``mirror()`` etc.
        Epochs are Julian years, e.g. ``2031.81``.

    Examples
    --------
    >>> from scopesim_templates.stellar import s301_orbit
    >>> orbit = s301_orbit(2)
    >>> t_peri = orbit.next_pericentre("2026-01-01")   # ~2031.81
    >>> rv = orbit.apparent_radial_velocity(t_peri)     # km/s
    """
    try:
        return S301_SOLUTIONS[int(solution)]
    except (KeyError, ValueError, TypeError) as err:
        raise ValueError(f"solution must be 1 or 2, not {solution!r}") \
            from err


def s301_spectrum(velocity, m_k=S301_MK, filter_name="Paranal/NACO.Ks",
                  template=S301_TEMPLATE, wave_range=(0.6, 2.7)):
    """Spectrum of S301 shifted to line-of-sight ``velocity``.

    ``m_k`` is applied *after* the shift, so it is the observed
    brightness. The wide default ``wave_range`` keeps the K band covered
    after shifts of up to ~8 per cent of c.

    Parameters
    ----------
    velocity : astropy.units.Quantity or float
        Line-of-sight (apparent) velocity; bare floats are taken as km/s.
    m_k : float
        Observed magnitude in ``filter_name`` (Vega).
    filter_name : str
        Filter for the normalisation (SVO ID or spextra name).
    template : str
        spextra template name.
    wave_range : tuple of float
        [um] Wavelength range to keep from the (very large) template.

    Returns
    -------
    spextra.Spextrum
    """
    velocity = u.Quantity(velocity, u.km / u.s)
    spec = Spextrum(template).cut(wave_range[0] * u.um, wave_range[1] * u.um)
    # Quantity, not float: Spextrum.redshift reads bare floats as m/s
    spec = spec.redshift(vel=velocity)
    return spec.scale_to_magnitude(m_k * u.mag, filter_curve=filter_name)


def sgr_a_star_spectrum(m_k=SGR_A_STAR_MK, alpha=SGR_A_STAR_ALPHA,
                        filter_name="Paranal/NACO.Ks",
                        wave_range=(0.6, 2.7)):
    """Featureless power-law continuum (F_nu ~ nu^alpha) for Sgr A*.

    At S301's pericentre the star is ~1.4 mas from Sgr A*, far inside the
    ELT's K-band diffraction limit, so the two blend into one point source.
    """
    waves = np.linspace(wave_range[0], wave_range[1], 2001) * u.um
    spec = Spextrum.powerlaw(alpha=alpha, x_0=2.2 * u.um,
                             amplitude=0 * u.ABmag, waves=waves,
                             filter_curve=filter_name)
    return spec.scale_to_magnitude(m_k * u.mag, filter_curve=filter_name)


@add_function_call_str
def s301(time,
         solution=1,
         with_sgr_a_star=True,
         m_k=S301_MK,
         sgr_a_star_mk=SGR_A_STAR_MK,
         relativistic=True,
         on_axis=False,
         filter_name="Paranal/NACO.Ks",
         template=S301_TEMPLATE):
    """
    Source of S301 (and optionally Sgr A*) at a given epoch.

    The star is placed at its projected offset from Sgr A* on the chosen
    orbit solution (Mang et al. 2026), and its spectrum is shifted to the
    apparent line-of-sight velocity at ``time`` (including transverse
    Doppler and gravitational redshift unless ``relativistic=False``).

    Parameters
    ----------
    time : astropy.time.Time, datetime, str, or float MJD
        Epoch. Same conventions as
        :func:`~scopesim_templates.stellar.galactic_centre`: bare numbers
        are MJD; for a decimal year use ``Time(2031.81, format="jyear")``
        or ``"J2031.81"``.
    solution : {1, 2}
        Which of the two mirror orbit solutions to use. They give the same
        sky position but opposite-signed (Newtonian) radial velocities.
    with_sgr_a_star : bool
        Add Sgr A*'s NIR counterpart as a featureless power law.
    m_k : float
        Observed K magnitude of S301 (19.3 +/- 0.3).
    sgr_a_star_mk : float
        Observed K magnitude of Sgr A*. 17.2 is a median quiescent state;
        flares reach ~15.
    relativistic : bool
        Include transverse-Doppler and gravitational redshift.
    on_axis : bool
        If False (default), coordinates are offsets from Sgr A*, so this
        Source can be added to
        :func:`~scopesim_templates.stellar.galactic_centre`. If True, S301
        is placed at (0, 0) and Sgr A* at its offset from S301 -- what an
        observer centring a slit on the star would do.
    filter_name : str
        Filter in which ``m_k`` and ``sgr_a_star_mk`` are defined.
    template : str
        spextra template for S301. Default ``"bosz/hr/f0v"``.

    Returns
    -------
    scopesim.Source
        Table columns ``x``, ``y`` [arcsec], ``ref``, ``weight``, ``rv``
        [km/s, apparent velocity applied to the spectrum] and ``name``.

    See Also
    --------
    s301_orbit : the orbit solutions themselves.
    galactic_centre : the Gillessen+2017 S-star field.

    Examples
    --------
    >>> from astropy.time import Time
    >>> from scopesim_templates.stellar import s301, s301_orbit
    >>> t = s301_orbit(1).next_pericentre("2026-01-01") + 1 / 365.25
    >>> src = s301(Time(t, format="jyear"), solution=1)
    """
    orbit = s301_orbit(solution)
    jyear = to_jyear(time)

    east, north = (float(np.ravel(v)[0]) for v in orbit.sky_position(jyear))
    rv = float(np.ravel(orbit.apparent_radial_velocity(
        jyear, relativistic=relativistic))[0])

    # mas -> arcsec; x = east (increasing RA), y = north, no RA flip
    star = np.array([east, north]) * 1e-3
    star_xy, hole_xy = ((np.zeros(2), -star) if on_axis
                        else (star, np.zeros(2)))

    spectra = [s301_spectrum(rv * u.km / u.s, m_k=m_k,
                             filter_name=filter_name, template=template)]
    positions, names, rvs = [star_xy], ["S301"], [rv]
    if with_sgr_a_star:
        spectra.append(sgr_a_star_spectrum(m_k=sgr_a_star_mk,
                                           filter_name=filter_name))
        positions.append(hole_xy)
        names.append("Sgr A*")
        rvs.append(0.0)

    x, y = np.transpose(positions)
    tbl = Table(
        names=["x", "y", "ref", "weight", "rv", "name"],
        data=[x, y, np.arange(len(names)), np.ones(len(names)), rvs, names],
        units=[u.arcsec, u.arcsec, None, None, u.km / u.s, None],
    )
    src = Source(spectra=spectra, table=tbl)
    src.meta.update({
        "description": f"S301 ({orbit.label}, Mang et al. 2026)",
        "epoch": str(time),
        "epoch_jyear": jyear,
        "solution": int(solution),
        "ra": RA_SGR_A_STAR,
        "dec": DEC_SGR_A_STAR,
        "on_axis": bool(on_axis),
        "rv_s301": rv,
        "relativistic": bool(relativistic),
        "sep_mas": float(np.hypot(east, north)),
        "m_bh": M_BH_SGR_A,
        "r0_pc": R0_SGR_A,
        "filter_name": filter_name,
    })
    return src
