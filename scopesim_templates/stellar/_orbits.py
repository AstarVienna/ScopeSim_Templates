# -*- coding: utf-8 -*-
"""Keplerian orbit machinery shared by the Galactic Centre and omega Cen templates.

Pure numpy; nothing here touches the network or any random number generator.

Contents
--------
* :func:`to_jyear` -- the epoch conventions used by all time-evolved
  templates (same as :func:`~scopesim_templates.stellar.galactic_centre`).
* :func:`solve_kepler` -- Halley solver for ``E - e sin E = M``, robust to
  the e = 0.98 orbits of S301.
* :class:`KeplerOrbit` -- one bound orbit about Sgr A* in the Galactic
  Centre literature convention (Gillessen et al. 2009, appendix), with
  sky positions, Newtonian and relativistic (apparent) radial velocities.
* :func:`sky_state` -- the vectorised core of :class:`KeplerOrbit`, used
  for many-star populations.
* :func:`kepler_propagate` -- universal-variable propagation of 3D
  position/velocity states (bound or unbound), used to time-evolve the
  omega Cen model.

Ported from the MAVIS 2026 workshop feasibility studies (``s301.py`` and
``omega_cen.py``).

Sgr A* parameters
-----------------
``M_BH_SGR_A`` and ``R0_SGR_A`` are from GRAVITY Collaboration 2022
(A&A 657, L12). Note that the S-star catalogue template
(:mod:`scopesim_templates.stellar._gillessen`) uses the older Gillessen et
al. 2017 values for consistency with its catalogue; the 0.5 per cent
difference in R0 is immaterial for imaging simulations.
"""

import warnings
from dataclasses import dataclass, replace
from datetime import datetime
from numbers import Real

import numpy as np
from astropy import units as u
from astropy.constants import G, c
from astropy.time import Time


__all__ = ["M_BH_SGR_A", "R0_SGR_A", "RA_SGR_A_STAR", "DEC_SGR_A_STAR",
           "KeplerOrbit", "to_jyear",
           "solve_kepler", "sky_state", "kepler_propagate"]


#: Mass of Sgr A* [Msun] (GRAVITY Collaboration 2022)
M_BH_SGR_A = 4.297e6
#: Distance to Sgr A* [pc] (GRAVITY Collaboration 2022)
R0_SGR_A = 8277.0
#: Position of Sgr A* (J2000) [deg]
RA_SGR_A_STAR = 266.41684
DEC_SGR_A_STAR = -29.00781

#: 1 au/yr in km/s
KMS_PER_AU_YR = (1 * u.au / u.yr).to_value(u.km / u.s)
#: speed of light in km/s
C_KMS = c.to_value(u.km / u.s)
#: G in pc^3 / Msun / yr^2 and in au^3 / Msun / yr^2
G_PCYR = G.to_value(u.pc**3 / u.Msun / u.yr**2)
G_AUYR = G.to_value(u.au**3 / u.Msun / u.yr**2)
#: Schwarzschild radius of 1 Msun, in au
R_S_PER_MSUN_AU = (2 * G * u.Msun / c**2).to_value(u.au)


# ---------------------------------------------------------------------------
# epochs
# ---------------------------------------------------------------------------

def to_jyear(time):
    """Convert an epoch to a Julian year (e.g. 2031.81).

    Accepts the same inputs as
    :func:`~scopesim_templates.stellar.galactic_centre`:
    ``astropy.time.Time``, ``datetime``, an ISO string (``"2031-10-22"``)
    or a Julian-epoch string (``"J2031.81"``), or a bare number, which is
    read as an **MJD**.

    A bare number below 10000 (i.e. an MJD before 1886) almost certainly
    means the caller passed a decimal year; that is still read as an MJD,
    as documented, but a warning is raised.
    """
    if isinstance(time, Time):
        return float(time.jyear)
    if isinstance(time, datetime):
        return float(Time(time).jyear)
    if isinstance(time, Real) and not isinstance(time, bool):
        if 0 < float(time) < 1e4:
            warnings.warn(
                f"Epoch {time!r} is interpreted as an MJD (year "
                f"{Time(float(time), format='mjd').jyear:.1f}). For a decimal "
                f"year pass e.g. Time({time!r}, format='jyear') or "
                f"'J{time}'.", UserWarning, stacklevel=3)
        return float(Time(float(time), format="mjd").jyear)
    return float(Time(time).jyear)


# ---------------------------------------------------------------------------
# Kepler's equation and the projected orbit
# ---------------------------------------------------------------------------

def solve_kepler(mean_anomaly, ecc, tol=1e-12, maxiter=60):
    """Solve ``E - e sin E = M`` for the eccentric anomaly, vectorised.

    Uses the Danby (1988) starting guess with Halley's method, which handles
    e = 0.983 without the stalling that plain Newton shows near pericentre.
    ``M`` is wrapped to ``[-pi, pi)``.

    Raises
    ------
    RuntimeError
        If the iteration has not converged after ``maxiter`` steps, rather
        than returning a silently wrong answer.
    """
    ecc = np.asarray(ecc, dtype=float)
    mean = np.mod(np.asarray(mean_anomaly, dtype=float) + np.pi,
                  2 * np.pi) - np.pi
    mean, ecc = np.broadcast_arrays(mean, ecc)
    ecc_anom = mean + 0.85 * ecc * np.sign(np.sin(mean))
    for _ in range(maxiter):
        sin_e, cos_e = np.sin(ecc_anom), np.cos(ecc_anom)
        f = ecc_anom - ecc * sin_e - mean
        if np.all(np.abs(f) < tol):
            break
        f1 = 1 - ecc * cos_e
        f2 = ecc * sin_e
        # Halley: one order better than Newton for the same derivatives
        ecc_anom = ecc_anom - 2 * f * f1 / (2 * f1**2 - f * f2)
    else:
        raise RuntimeError("Kepler solver did not converge "
                           f"(max e = {np.max(ecc)})")
    return ecc_anom


def sky_state(a_mas, ecc, incl, node, argp, t_peri, period, epoch,
              distance_pc=R0_SGR_A, central_mass=M_BH_SGR_A):
    """Projected position and line-of-sight motion of Keplerian orbits.

    All orbital elements broadcast against each other (one value per star)
    and against ``epoch``. Angles are in degrees, in the Galactic Centre
    convention: ``node`` is the position angle of the ascending node
    measured North through East, the ascending node is crossed moving
    *away* from the observer, and positive radial velocity is receding.

    Parameters
    ----------
    a_mas : array_like
        Semi-major axis [mas].
    ecc, incl, node, argp : array_like
        Eccentricity and the three angles [deg].
    t_peri, period : array_like
        Epoch of pericentre and period [Julian years].
    epoch : float or array_like
        Julian year(s) at which to evaluate the orbits.
    distance_pc, central_mass : float
        Distance [pc] and central mass [Msun], used only for the speed and
        the relativistic terms.

    Returns
    -------
    dict
        ``east``, ``north`` [mas], ``radius`` [au], ``v_r`` (Newtonian
        line-of-sight velocity), ``speed`` (vis-viva) [km/s].
    """
    a_mas, ecc, incl, node, argp, t_peri, period, epoch = np.broadcast_arrays(
        *(np.asarray(v, dtype=float) for v in
          (a_mas, ecc, incl, node, argp, t_peri, period, epoch)))
    mean = 2 * np.pi * (epoch - t_peri) / period
    ecc_anom = solve_kepler(mean, ecc)
    cos_e, sin_e = np.cos(ecc_anom), np.sin(ecc_anom)
    nu = np.arctan2(np.sqrt(1 - ecc**2) * sin_e, cos_e - ecc)
    rel_radius = 1 - ecc * cos_e                     # r / a

    theta = nu + np.radians(argp)
    om, inc = np.radians(node), np.radians(incl)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    scale = rel_radius * a_mas
    north = scale * (np.cos(om) * cos_t - np.sin(om) * sin_t * np.cos(inc))
    east = scale * (np.sin(om) * cos_t + np.cos(om) * sin_t * np.cos(inc))

    # 1 arcsec at 1 pc is 1 au, by construction of the parsec
    a_au = a_mas * 1e-3 * distance_pc
    k_amp = (2 * np.pi * a_au / period * np.sin(inc)
             / np.sqrt(1 - ecc**2)) * KMS_PER_AU_YR
    v_r = k_amp * (np.cos(theta) + ecc * np.cos(np.radians(argp)))

    radius = a_au * rel_radius
    speed = np.sqrt(G_AUYR * central_mass * (2 / radius - 1 / a_au))
    return {"east": east, "north": north, "radius": radius,
            "v_r": v_r, "speed": speed * KMS_PER_AU_YR}


def apparent_velocity(v_r, speed, radius_au, central_mass=M_BH_SGR_A):
    """Velocity a spectroscopist would report, ``c z``, in km/s.

    Leading-order combination of Doppler, transverse-Doppler and
    gravitational redshift::

        1 + z = (1 + v_r/c) / sqrt(1 - v^2/c^2) / sqrt(1 - R_s/r)

    For S301 at pericentre the two relativistic terms add ~2200 km/s.
    """
    beta_r = np.asarray(v_r) / C_KMS
    beta = np.asarray(speed) / C_KMS
    r_ratio = R_S_PER_MSUN_AU * central_mass / np.asarray(radius_au)
    one_plus_z = (1 + beta_r) / np.sqrt(1 - beta**2) / np.sqrt(1 - r_ratio)
    return (one_plus_z - 1) * C_KMS


@dataclass(frozen=True)
class KeplerOrbit:
    """A bound orbit about Sgr A* (or any point mass).

    Angles in degrees, ``t_peri`` and ``period`` in Julian years, ``a`` in
    milliarcseconds on the sky. Conventions as in :func:`sky_state`.
    Epoch arguments of the methods are Julian years (floats); use
    :func:`to_jyear` to convert ``Time`` objects.
    """

    label: str
    a: float          # mas
    e: float
    incl: float       # deg
    node: float       # deg, Omega
    argp: float       # deg, omega
    t_peri: float     # Julian year
    period: float     # yr
    distance_pc: float = R0_SGR_A
    central_mass: float = M_BH_SGR_A

    @property
    def a_au(self):
        """Semi-major axis [au]."""
        return self.a * 1e-3 * self.distance_pc

    @property
    def r_peri_au(self):
        """Pericentre distance [au]."""
        return self.a_au * (1 - self.e)

    @property
    def r_peri_rs(self):
        """Pericentre distance in Schwarzschild radii."""
        return self.r_peri_au / (R_S_PER_MSUN_AU * self.central_mass)

    @property
    def v_peri(self):
        """Orbital speed at pericentre [km/s] (Newtonian, vis-viva)."""
        return float(np.sqrt(G_AUYR * self.central_mass
                             * (2 / self.r_peri_au - 1 / self.a_au))
                     * KMS_PER_AU_YR)

    @property
    def kepler_period(self):
        """Period [yr] implied by ``a`` and ``central_mass`` (Kepler III)."""
        return float(np.sqrt(self.a_au**3 / self.central_mass))

    @property
    def precession_per_orbit(self):
        """Schwarzschild periapsis advance per orbit [deg] (not applied)."""
        r_s = R_S_PER_MSUN_AU * self.central_mass
        return float(np.degrees(3 * np.pi * r_s
                                / (self.a_au * (1 - self.e**2))))

    def pericentre(self, n=0):
        """Julian year of the ``n``-th pericentre after ``t_peri``."""
        return self.t_peri + n * self.period

    def next_pericentre(self, after):
        """Julian year of the first pericentre after the epoch ``after``."""
        after = after if isinstance(after, Real) else to_jyear(after)
        return self.pericentre(int(np.ceil((after - self.t_peri)
                                           / self.period)))

    def mirror(self):
        """The mirror solution: same sky ellipse, opposite radial velocity."""
        return replace(self, label=self.label + " (mirrored)",
                       node=(self.node + 180) % 360,
                       argp=(self.argp + 180) % 360)

    def _state(self, epoch):
        return sky_state(self.a, self.e, self.incl, self.node, self.argp,
                         self.t_peri, self.period, epoch,
                         self.distance_pc, self.central_mass)

    def sky_position(self, epoch):
        """``(east, north)`` offsets from the central mass [mas]."""
        state = self._state(epoch)
        return state["east"], state["north"]

    def separation(self, epoch):
        """Projected separation from the central mass [mas]."""
        return np.hypot(*self.sky_position(epoch))

    def radial_velocity(self, epoch):
        """Newtonian line-of-sight velocity [km/s], positive = receding."""
        return self._state(epoch)["v_r"]

    def apparent_radial_velocity(self, epoch, relativistic=True):
        """Observed ``c z`` [km/s], see :func:`apparent_velocity`."""
        state = self._state(epoch)
        if not relativistic:
            return state["v_r"]
        return apparent_velocity(state["v_r"], state["speed"],
                                 state["radius"], self.central_mass)


# ---------------------------------------------------------------------------
# universal-variable propagation of 3D states
# ---------------------------------------------------------------------------

def _stumpff(z):
    """Stumpff functions c2(z) and c3(z), stable across the sign of z."""
    c2, c3 = np.empty_like(z), np.empty_like(z)
    pos, neg = z > 1e-8, z < -1e-8
    small = ~(pos | neg)
    root = np.sqrt(np.abs(z))

    c2[pos] = (1 - np.cos(root[pos])) / z[pos]
    c3[pos] = (root[pos] - np.sin(root[pos])) / root[pos] ** 3
    c2[neg] = (np.cosh(root[neg]) - 1) / -z[neg]
    c3[neg] = (np.sinh(root[neg]) - root[neg]) / root[neg] ** 3
    c2[small] = 0.5 - z[small] / 24
    c3[small] = 1 / 6 - z[small] / 120
    return c2, c3


def kepler_propagate(pos, vel, mu, dt, tol=1e-11, maxiter=100):
    """Advance positions and velocities on Keplerian orbits by ``dt``.

    Universal-variable formulation with Laguerre iteration, so elliptic,
    parabolic and hyperbolic orbits go through the same code path. Units
    are arbitrary but consistent, e.g. ``pos`` in pc, ``vel`` in pc/yr,
    ``mu = G M`` in pc^3/yr^2 and ``dt`` in yr.

    Parameters
    ----------
    pos, vel : ndarray, shape (n, 3)
    mu : float or ndarray, shape (n,)
        Gravitational parameter of the central mass seen by each body.
    dt : float

    Returns
    -------
    new_pos, new_vel : ndarray, shape (n, 3)
        New arrays; the inputs are not modified.
    """
    pos = np.asarray(pos, dtype=float)
    vel = np.asarray(vel, dtype=float)
    if not len(pos) or not dt:
        return pos.copy(), vel.copy()
    mu = np.broadcast_to(np.asarray(mu, dtype=float), (len(pos),))

    r0 = np.linalg.norm(pos, axis=1)
    v0sq = np.einsum("ij,ij->i", vel, vel)
    sqrt_mu = np.sqrt(mu)
    sigma0 = np.einsum("ij,ij->i", pos, vel) / sqrt_mu
    alpha = 2.0 / r0 - v0sq / mu          # = 1/a

    # For bound orbits, fold dt into a single radial period: keeps the
    # initial guess good no matter how many revolutions have elapsed.
    dt_eff = np.full_like(r0, float(dt))
    bound = alpha > 1e-14
    period = np.where(bound,
                      2 * np.pi / (sqrt_mu * np.abs(alpha) ** 1.5), np.inf)
    dt_eff[bound] = np.remainder(dt, period[bound])

    chi = np.where(bound, sqrt_mu * dt_eff * alpha, sqrt_mu * dt_eff / r0)
    target = sqrt_mu * dt_eff

    for _ in range(maxiter):
        z = alpha * chi**2
        c2, c3 = _stumpff(z)
        f = (r0 * chi + sigma0 * chi**2 * c2
             + (1 - alpha * r0) * chi**3 * c3 - target)
        df = (r0 + sigma0 * chi * (1 - z * c3)
              + (1 - alpha * r0) * chi**2 * c2)
        d2f = sigma0 * (1 - z * c2) + (1 - alpha * r0) * chi * (1 - z * c3)

        # Laguerre's method (n = 5): far more robust than Newton here
        radicand = np.abs(16 * df**2 - 20 * f * d2f)
        sign = np.where(df >= 0.0, 1.0, -1.0)
        step = 5 * f / (df + sign * np.sqrt(radicand))
        chi = chi - step
        if np.all(np.abs(step) <= tol * (np.abs(chi) + 1.0)):
            break
    else:
        warnings.warn("Kepler solver did not fully converge; worst residual "
                      f"{np.abs(f).max():.3e}", stacklevel=3)

    z = alpha * chi**2
    c2, c3 = _stumpff(z)
    lf = 1 - chi**2 / r0 * c2
    lg = dt_eff - chi**3 / sqrt_mu * c3
    new_pos = lf[:, None] * pos + lg[:, None] * vel

    r1 = np.linalg.norm(new_pos, axis=1)
    lfdot = sqrt_mu / (r0 * r1) * chi * (z * c3 - 1)
    lgdot = 1 - chi**2 / r1 * c2
    new_vel = lfdot[:, None] * pos + lgdot[:, None] * vel
    return new_pos, new_vel
