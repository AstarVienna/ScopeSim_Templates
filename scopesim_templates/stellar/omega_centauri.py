# -*- coding: utf-8 -*-
"""Toy model of omega Centauri (NGC 5139) with a central intermediate-mass
black hole, and Keplerian time evolution of the stellar positions.

This is the detailed, literature-calibrated counterpart of the generic
:func:`~scopesim_templates.stellar.globular_cluster` template (Salpeter IMF,
random Kepler orbits with log-uniform semi-major axes, fixed star count).
Use :func:`omega_cen` when the structure, luminosity function and
kinematics should look like the real cluster.

The cluster is built in two independent halves:

**The population** -- a Kroupa (2001) IMF is sampled down to 0.01 Msun and
truncated at the main-sequence turn-off of an old, metal-poor population.
Stars above the turn-off are dark remnants and are dropped; an evolved
(RGB) component is added separately, since those stars dominate the light of
a globular cluster. Spectra come from the Pickles library (which covers the
optical *and* the near-IR, so one ``Source`` serves MAVIS and MICADO alike),
preferring the weak-lined metal-poor variants.

**The dynamics** -- stars are drawn from a spherical, isotropic "analytic
King" (modified Hubble) density profile, given velocities from the isotropic
Jeans equation *including* the IMBH, and propagated on Keplerian orbits about
the central mass by :func:`evolve_omega_cen`.

Structural parameters are from Harris (1996, 2010 ed.), Baumgardt & Hilker
(2018) and Baumgardt & Vasiliev (2021); the IMBH is the one inferred by
Haeberle et al. (2024, Nature 631, 285), whose *lower limit* is 8200 Msun.

How closely it matches the real thing
-------------------------------------
Checked in the workshop against the literature at the default settings:

===========================  ==========  ============================
quantity                     model       observed
===========================  ==========  ============================
central density log rho_0    3.27        3.15-3.3 (Harris)
central mu_V [mag/arcsec2]   17.1        16.77 (Harris)
main-sequence turn-off       V = 17.5    V = 17.5 (by construction)
RGB tip                      V = 11.2    V ~ 11.5
stars brighter than M_V = 0  ~2300       ~1500-2500
sigma_1D at 1 pc [km/s]      20.8        17-22
IMBH influence radius        2.3 pc      (89 arcsec at 5.43 kpc)
===========================  ==========  ============================

Caveats
-------
* No mass segregation, no binaries, no rotation, no anisotropy.
* The density profile has a flat core. An IMBH is expected to grow a
  Bahcall & Wolf (1976) rho ~ r^-7/4 cusp; omega Cen shows no such cusp
  observationally, so ``cusp_slope`` defaults to 0 but is available.
* Star counts use the present-day dynamical mass with an *initial* mass
  function, which over-counts stars by a few tens of per cent.
* There is no sub-giant or horizontal branch, and the main sequence comes
  from a zero-age table (corrected near the turn-off via a "photometric
  mass", see :func:`_photometric_mass`), so it supplies only ~20 per cent of
  the V light instead of ~55 per cent. The evolved component, normalised to
  the total mass-to-light ratio, absorbs the difference: the bright end of
  the luminosity function is right, but stars within a magnitude of the
  turn-off are labelled "RGB" more often than they should be.
* Velocities near the IMBH are truncated at the local escape speed: a cored
  profile is not in equilibrium around a point mass.
* Spectra are shared per template (not RV-shifted per star); the 3D
  velocity is in the ``vx``, ``vy``, ``vz`` columns.

Differences from the workshop code (``omega_cen.py``)
-----------------------------------------------------
* The IMF is sampled with an exact inverse-CDF sampler on an
  ``np.random.Generator`` instead of ``scopesim_targets``'
  ``BrokenPowerlaw`` + ``scipy.stats.sampling``; absolute magnitudes come
  from the package's own Mamajek table instead of ``StellarParameters``.
* Dwarfs cooler than the coolest Pickles dwarf (M6V, ~0.1 Msun) use the M6V
  template instead of a blackbody. They are ~20 mag below the turn-off.
* Epochs follow the package convention (``Time``, ISO / ``"J2030"``
  strings, or MJD), not bare decimal years.
"""

import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from astropy import units as u
from astropy.constants import G as _G
from astropy.table import Table

import pyckles
from spextra import Spextrum

from ..rc import Source
from ..utils.general_utils import add_function_call_str
from ._orbits import kepler_propagate, to_jyear
from ._population import (sample_broken_powerlaw, broken_powerlaw_mean,
                          dwarf_absmag, dwarf_table, seed_for_meta)


__all__ = ["omega_cen", "evolve_omega_cen"]


# ---------------------------------------------------------------------------
# omega Cen, and other constants
# ---------------------------------------------------------------------------

RA_OMEGA_CEN = 201.69683    # [deg] J2000, Harris (13:26:47.24)
DEC_OMEGA_CEN = -47.47958   # [deg] J2000, Harris (-47:28:46.5)

DISTANCE = 5430.0        # [pc] Baumgardt & Vasiliev 2021, cf. oMEGACat VI
CLUSTER_MASS = 3.55e6    # [Msun] Baumgardt & Hilker 2018
R_CORE = 3.74            # [pc] Harris r_c = 2.37 *arcmin* at 5.43 kpc
R_TIDE = 90.0            # [pc] Harris r_t = 57.0 arcmin at 5.43 kpc
IMBH_MASS = 1e5          # [Msun] cf. Haeberle et al. 2024 (>= 8.2e3)

TURNOFF_MASS = 0.83      # [Msun] main-sequence turn-off, 12 Gyr, [Fe/H]=-1.6
TURNOFF_ABSMAG = 3.85    # [mag] observed M_V of the turn-off (V ~ 17.5)
UNEVOLVED_MASS = 0.5     # [Msun] below this, 12 Gyr stars are still ZAMS
MIN_MASS = 0.01          # [Msun] below the deuterium-burning limit
IMF_TOP = 60.0           # [Msun]

#: [Msun/Lsun] omega Cen in V (Baumgardt & Hilker 2018). Fixes how many
#: evolved stars there are: they supply the light the main sequence lacks.
MASS_TO_LIGHT = 2.6

#: All magnitude bookkeeping is in Johnson V; the templates are real
#: spectra, so fluxes in every other band follow from them.
REF_FILTER = "Generic/Johnson.V"
SOLAR_ABSMAG_V = 4.83

#: Luminosity function of the evolved component, dN/dM_V ~ 10^(slope*M_V)
RGB_LF_SLOPE = 0.35
RGB_MAG_RANGE = (-2.5, 3.5)   # [mag] M_V from the tip down to sub-giants

#: Metal-poor RGB: absolute V magnitude -> Pickles spectrum. Metal-poor
#: giants are hotter than solar-metallicity ones at fixed luminosity, hence
#: the early types; weak-lined ("_W") variants are used where they exist.
RGB_SEQUENCE = [
    (3.5, "G8III_W"), (2.5, "K0III_W"), (1.5, "K1III_W"), (0.5, "K2III_W"),
    (-0.5, "K3III_W"), (-1.0, "K4III_W"), (-1.5, "K5III"), (-2.0, "M0III"),
    (-2.5, "M2III"),
]

DEFAULT_EPOCH = "J2030.0"

# G in the units used below: pc, Msun, km/s and pc/yr
G_KMS = _G.to_value(u.pc * u.km**2 / u.s**2 / u.Msun)
G_PCYR = _G.to_value(u.pc**3 / u.Msun / u.yr**2)
KMS_PER_PCYR = (1 * u.pc / u.yr).to_value(u.km / u.s)
RAD2ARCSEC = u.rad.to(u.arcsec)

#: columns that :func:`evolve_omega_cen` needs
_STATE_COLUMNS = ("x_pc", "y_pc", "z_pc", "vx", "vy", "vz")


# ---------------------------------------------------------------------------
# structure: density, enclosed mass, velocity dispersion
# ---------------------------------------------------------------------------

def _density_profile(r, r_core, cusp_slope=0.0, r_cusp=None):
    """Unnormalised 3D density: analytic King, optionally with a cusp."""
    rho = (1.0 + (r / r_core) ** 2) ** -1.5
    if cusp_slope:
        rho = rho * np.maximum(1.0, (r / r_cusp) ** -cusp_slope)
    return rho


class _ClusterPotential:
    """Enclosed mass, potential and Jeans dispersion of the toy cluster.

    Everything is tabulated on a log-spaced radial grid so that the density
    profile can be changed (e.g. given a cusp) without analytic rework.
    """

    def __init__(self, mass=CLUSTER_MASS, r_core=R_CORE, r_tide=R_TIDE,
                 imbh_mass=IMBH_MASS, cusp_slope=0.0, n_grid=6001):
        self.mass = float(mass)
        self.r_core = float(r_core)
        self.r_tide = float(r_tide)
        self.imbh_mass = float(imbh_mass)
        self.cusp_slope = float(cusp_slope)

        # the cusp, if any, is anchored at the IMBH radius of influence
        r_cusp = self._influence_radius_flat_core()

        self.r = np.geomspace(1e-6 * self.r_core, self.r_tide, n_grid)
        rho = _density_profile(self.r, self.r_core, self.cusp_slope, r_cusp)

        # enclosed cluster mass, normalised so that M(r_tide) = mass
        m_cum = cumulative_trapezoid(4 * np.pi * self.r**2 * rho, self.r,
                                     initial=0.0)
        self.rho = rho * self.mass / m_cum[-1]
        self.m_encl = m_cum * self.mass / m_cum[-1]
        self.r_cusp = r_cusp
        self._setup_kinematics()

    def _influence_radius_flat_core(self):
        """Radius where the flat-core cluster mass equals the IMBH mass."""
        x_t = self.r_tide / self.r_core
        rho0 = self.mass / (4 * np.pi * self.r_core**3 * (
            np.arcsinh(x_t) - x_t / np.hypot(1.0, x_t)))
        return (3 * self.imbh_mass / (4 * np.pi * rho0)) ** (1 / 3)

    def _setup_kinematics(self):
        m_tot = self.imbh_mass + self.m_encl

        # isotropic Jeans equation, integrated inwards from sigma(r_tide)=0
        integrand = self.rho * G_KMS * m_tot / self.r**2
        outward = cumulative_trapezoid(integrand, self.r, initial=0.0)
        self.sigma = np.sqrt(np.maximum(outward[-1] - outward, 0.0)
                             / self.rho)

        # potential, for the escape speed: -G[M(r)/r + 4 pi int rho s ds]
        shells = cumulative_trapezoid(4 * np.pi * self.rho * self.r, self.r,
                                      initial=0.0)
        self.phi = -G_KMS * (m_tot / self.r + shells[-1] - shells)

        # cumulative number fraction, for inverse-transform sampling
        self.mass_frac = self.m_encl / self.m_encl[-1]

    def enclosed_mass(self, r):
        """Cluster mass (no IMBH) inside radius ``r`` [pc], in Msun."""
        return np.interp(r, self.r, self.m_encl)

    def central_mass(self, r):
        """IMBH plus cluster mass inside ``r`` [pc], in Msun."""
        return self.imbh_mass + self.enclosed_mass(r)

    def dispersion(self, r):
        """Isotropic 1D velocity dispersion at ``r`` [pc], in km/s."""
        return np.interp(r, self.r, self.sigma)

    def escape_speed(self, r):
        """Local escape speed at ``r`` [pc], in km/s."""
        return np.sqrt(-2 * np.interp(r, self.r, self.phi))

    def sample_radii(self, n, rng):
        """Draw ``n`` 3D radii [pc] from the density profile."""
        return np.interp(rng.random(n), self.mass_frac, self.r)


# ---------------------------------------------------------------------------
# spatial and velocity sampling
# ---------------------------------------------------------------------------

def _sample_cylinder(potential, rng, half_size, n_draw=None, n_keep=None,
                     batch=1_000_000, max_batches=1000):
    """Sample 3D positions [pc] whose projection lands in the square field.

    Exactly one of ``n_draw`` (draw this many cluster stars and keep those
    in the field) or ``n_keep`` (draw until this many are in the field) must
    be given. Stars are placed isotropically, so the line-of-sight extent is
    the full cluster.
    """
    kept, n_drawn, n_kept = [], 0, 0
    for _ in range(max_batches):
        if n_draw is not None:
            size = int(min(batch, n_draw - n_drawn))
        else:
            # over-draw a little to avoid a long tail of tiny batches
            rate = n_kept / n_drawn if n_drawn else 0.0
            size = int(min(batch, max(1000, 1.2 * (n_keep - n_kept)
                                      / max(rate, 1e-6))))
        if size <= 0:
            break

        radius = potential.sample_radii(size, rng)
        costheta = rng.uniform(-1.0, 1.0, size)
        sintheta = np.sqrt(1.0 - costheta**2)
        phi = rng.uniform(0.0, 2 * np.pi, size)

        xyz = np.stack([radius * sintheta * np.cos(phi),
                        radius * sintheta * np.sin(phi),
                        radius * costheta], axis=1)
        hit = ((np.abs(xyz[:, 0]) < half_size)
               & (np.abs(xyz[:, 1]) < half_size))

        n_drawn += size
        kept.append(xyz[hit])
        n_kept += int(hit.sum())

        if n_draw is not None and n_drawn >= n_draw:
            break
        if n_keep is not None and n_kept >= n_keep:
            break
    else:
        if n_keep is not None and n_kept < n_keep:
            warnings.warn(f"Only placed {n_kept} of {n_keep} stars in the "
                          "field before giving up; the field is very small "
                          "compared with the cluster.", stacklevel=3)

    xyz = np.concatenate(kept) if kept else np.empty((0, 3))
    if n_keep is not None:
        xyz = xyz[:n_keep]
    return xyz, n_drawn


def _sample_velocities(potential, xyz, rng, maxiter=100):
    """Isotropic Gaussian velocities [km/s] from the Jeans dispersion.

    Speeds above the local escape speed are redrawn: around a point mass
    the Jeans equation for a *cored* profile gives sigma_r = v_circ, so an
    untruncated Gaussian would unbind most of the innermost stars.
    """
    radius = np.linalg.norm(xyz, axis=1)
    sigma = potential.dispersion(radius)
    v_esc = potential.escape_speed(radius)

    vel = rng.normal(0.0, 1.0, xyz.shape) * sigma[:, None]
    unbound = np.linalg.norm(vel, axis=1) >= v_esc
    for _ in range(maxiter):
        if not unbound.any():
            break
        redraw = rng.normal(0.0, 1.0, (int(unbound.sum()), 3))
        vel[unbound] = redraw * sigma[unbound, None]
        unbound[unbound] = (np.linalg.norm(vel[unbound], axis=1)
                            >= v_esc[unbound])
    if unbound.any():  # stubborn cases: rescale onto a bound speed
        speed = np.linalg.norm(vel[unbound], axis=1)
        vel[unbound] *= (0.9 * v_esc[unbound] / speed)[:, None]
    return vel


# ---------------------------------------------------------------------------
# photometry
# ---------------------------------------------------------------------------

def _photometric_mass(masses, turnoff_mass, turnoff_absmag):
    """Map true masses onto the ZAMS masses with the same photometry.

    The Mamajek table describes zero-age, solar-metallicity dwarfs, so it
    puts a 0.83 Msun star at M_V ~ 6. In a 12 Gyr, [Fe/H] ~ -1.6 population
    that star sits at the turn-off, at M_V ~ 3.85: hotter, bluer and two
    magnitudes brighter. Looking stars up at an inflated "photometric mass"
    recovers luminosity *and* colour from the observed turn-off magnitude.
    The correction tapers to zero at ``UNEVOLVED_MASS``.
    """
    grid = np.linspace(0.08, 2.0, 400)
    absmag = dwarf_absmag(grid, "Mv")
    order = np.argsort(absmag)
    zams_mass = np.interp(turnoff_absmag, absmag[order], grid[order])

    ramp = np.clip((masses - UNEVOLVED_MASS)
                   / (turnoff_mass - UNEVOLVED_MASS), 0.0, 1.0)
    return masses * (1.0 + (zams_mass / turnoff_mass - 1.0) * ramp)


def _n_giants(mass, mass_to_light, ms_light):
    """How many evolved stars the cluster needs to reach its observed M/L.

    The main sequence supplies ``ms_light`` [Lsun, whole cluster] and the
    RGB makes up the rest of ``mass / mass_to_light``.
    """
    lo, hi = RGB_MAG_RANGE
    # mean luminosity of one evolved star, over dN/dM ~ 10^(slope * M)
    weighted = quad(lambda m: 10 ** (-0.4 * (m - SOLAR_ABSMAG_V)
                                     + RGB_LF_SLOPE * m), lo, hi)[0]
    norm = quad(lambda m: 10 ** (RGB_LF_SLOPE * m), lo, hi)[0]
    missing_light = mass / mass_to_light - ms_light
    return max(0.0, missing_light) / (weighted / norm)


def _sample_giants(n, rng, metal_poor=True):
    """Absolute V magnitudes and Pickles template names of ``n`` RGB stars."""
    slope, (mag_hi, mag_lo) = RGB_LF_SLOPE, RGB_MAG_RANGE
    lo, hi = 10 ** (slope * mag_hi), 10 ** (slope * mag_lo)
    absmag = np.log10(lo + rng.random(n) * (hi - lo)) / slope

    mags = np.array([mag for mag, _ in RGB_SEQUENCE])[::-1]
    names = [name for _, name in RGB_SEQUENCE][::-1]
    if not metal_poor:
        names = [name.removesuffix("_W") for name in names]
    index = np.rint(np.interp(absmag, mags,
                              np.arange(len(mags)))).astype(int)
    return absmag, np.array(names, dtype=object)[index]


def _pickles_dwarfs(pickles_names, metal_poor=True):
    """Map available Pickles dwarf templates to their Mamajek mass.

    Returns ``(names, masses)`` sorted by mass, one template per spectral
    type, preferring the weak-lined ``_W`` variant if ``metal_poor``.
    """
    mamajek = dwarf_table()
    spt_mass = {str(spt).strip(): float(m)
                for spt, m in zip(mamajek["SpT"], mamajek["Msun"])}
    variants = {}
    for name in map(str, pickles_names):
        base, _, variant = name.partition("_")
        if base.endswith("V") and not base.endswith("IV"):
            variants.setdefault(base, {})[variant] = name

    chosen = []
    for base, avail in variants.items():
        # Pickles writes M2.5V as "M25V"
        spt = base if len(base) <= 3 else f"{base[:2]}.{base[2:]}"
        if spt not in spt_mass:
            continue
        if metal_poor and "W" in avail:
            name = avail["W"]
        else:
            name = avail.get("", next(iter(avail.values())))
        chosen.append((name, spt_mass[spt]))
    names, masses = zip(*sorted(chosen, key=lambda nm: nm[1]))
    return list(names), np.array(masses)


def _assign_dwarf_templates(phot_mass, pickles_names, metal_poor=True):
    """Name of the Pickles dwarf closest in (log) mass to each star."""
    names, masses = _pickles_dwarfs(pickles_names, metal_poor)
    log_m = np.log10(np.clip(phot_mass, masses[0], masses[-1]))
    idx = np.abs(log_m[:, None] - np.log10(masses)[None, :]).argmin(axis=1)
    return np.array(names, dtype=object)[idx]


# ---------------------------------------------------------------------------
# the cluster
# ---------------------------------------------------------------------------

def _build_table(xyz, vel, mag, absmag, mass, phase, spec_types, distance):
    """Point-source table with positions, photometry and the 3D state."""
    to_mas_yr = 1.0 / (4.740470446 * distance / 1000.0)  # per km/s
    return Table({
        "x": xyz[:, 0] / distance * RAD2ARCSEC * u.arcsec,
        "y": xyz[:, 1] / distance * RAD2ARCSEC * u.arcsec,
        "mag": mag * u.mag,
        "absmag": absmag * u.mag,
        "mass": mass * u.Msun,
        "phase": np.asarray(phase, dtype="<U3"),
        "spec_types": np.asarray(spec_types, dtype=str),
        "x_pc": xyz[:, 0] * u.pc,
        "y_pc": xyz[:, 1] * u.pc,
        "z_pc": xyz[:, 2] * u.pc,
        "vx": vel[:, 0] * u.km / u.s,
        "vy": vel[:, 1] * u.km / u.s,
        "vz": vel[:, 2] * u.km / u.s,
        "pm_x": vel[:, 0] * to_mas_yr * u.mas / u.yr,
        "pm_y": vel[:, 1] * to_mas_yr * u.mas / u.yr,
    })


def _omega_cen_catalogue(fov=60.0, imbh_mass=IMBH_MASS, mass=CLUSTER_MASS,
                         distance=DISTANCE, r_core=R_CORE, r_tide=R_TIDE,
                         cusp_slope=0.0, turnoff_mass=TURNOFF_MASS,
                         turnoff_absmag=TURNOFF_ABSMAG, min_mass=MIN_MASS,
                         giants=True, mass_to_light=MASS_TO_LIGHT,
                         metal_poor=True, mag_limit=None,
                         max_stars=2_000_000, seed=None):
    """Sample the cluster stars in the field (no spectra, no network).

    Random draws happen in a fixed order -- positions, masses, giants, then
    velocities -- so that at fixed ``seed`` changing ``imbh_mass`` changes
    only the velocities. Main-sequence rows have an empty ``spec_types``,
    filled in later by :func:`omega_cen`.
    """
    rng = np.random.default_rng(seed)
    potential = _ClusterPotential(mass, r_core, r_tide, imbh_mass, cusp_slope)

    half_size = np.deg2rad(fov / 3600.0) / 2 * distance  # [pc]
    mean_mass = broken_powerlaw_mean(min_mass, IMF_TOP)
    n_cluster = int(round(mass / mean_mass))

    # 1. positions of every cluster star projected into the field
    xyz, n_drawn = _sample_cylinder(potential, rng, half_size,
                                    n_draw=n_cluster)
    if not len(xyz):
        raise ValueError(
            f"No cluster stars fall in the {fov} arcsec field. It is too "
            f"small for a {mass:.3g} Msun cluster with r_core = {r_core} pc.")
    if len(xyz) > max_stars:
        raise ValueError(
            f"{len(xyz)} stars fall in the {fov} arcsec field, more than "
            f"max_stars={max_stars}. Reduce fov or raise max_stars.")
    field_fraction = len(xyz) / n_drawn

    # 2. masses, keeping only stars still on the main sequence
    masses = sample_broken_powerlaw(len(xyz), min_mass, IMF_TOP, rng)
    alive = masses <= turnoff_mass
    xyz, masses = xyz[alive], masses[alive]

    phot_mass = _photometric_mass(masses, turnoff_mass, turnoff_absmag)
    absmag = dwarf_absmag(phot_mass, "Mv")
    phase = np.full(len(masses), "MS", dtype="<U3")
    spec_types = np.full(len(masses), "", dtype=object)

    # 3. the evolved stars, which carry most of the light. Their number
    #    comes from the light the main sequence fails to account for,
    #    scaled from the field up to the whole cluster.
    ms_light = np.sum(10 ** (-0.4 * (absmag - SOLAR_ABSMAG_V))) \
        / field_fraction
    n_giants_total = _n_giants(mass, mass_to_light, ms_light)
    n_giants = int(round(n_giants_total * field_fraction)) if giants else 0

    if n_giants:
        giant_xyz, _ = _sample_cylinder(potential, rng, half_size,
                                        n_keep=n_giants)
        giant_absmag, giant_names = _sample_giants(len(giant_xyz), rng,
                                                   metal_poor)
        xyz = np.concatenate([xyz, giant_xyz])
        absmag = np.concatenate([absmag, giant_absmag])
        masses = np.concatenate([masses, np.full(len(giant_xyz),
                                                 turnoff_mass)])
        phot_mass = np.concatenate([phot_mass, np.full(len(giant_xyz),
                                                       np.nan)])
        phase = np.concatenate([phase, np.full(len(giant_xyz), "RGB")])
        spec_types = np.concatenate([spec_types, giant_names])

    # 4. velocities, then assemble
    vel = _sample_velocities(potential, xyz, rng)
    distmod = 5 * np.log10(distance) - 5
    apparent = absmag + distmod

    if mag_limit is not None:
        bright = apparent <= mag_limit
        xyz, vel, apparent, absmag = (xyz[bright], vel[bright],
                                      apparent[bright], absmag[bright])
        masses, phot_mass = masses[bright], phot_mass[bright]
        phase, spec_types = phase[bright], spec_types[bright]

    tbl = _build_table(xyz, vel, apparent, absmag, masses, phase,
                       spec_types, distance)
    tbl["phot_mass"] = phot_mass * u.Msun
    tbl.meta.update({
        "imbh_mass": float(imbh_mass),
        "cluster_mass": float(mass),
        "distance_pc": float(distance),
        "r_core_pc": float(r_core),
        "r_tide_pc": float(r_tide),
        "cusp_slope": float(cusp_slope),
        "turnoff_mass": float(turnoff_mass),
        "turnoff_absmag": float(turnoff_absmag),
        "mass_to_light": float(mass_to_light),
        "n_giants_cluster": float(n_giants_total),
        "fov": float(fov),
        "n_cluster": n_cluster,
        "field_fraction": float(field_fraction),
        "sigma_1pc": float(potential.dispersion(1.0)),
        "r_influence_pc": float(potential.r_cusp),
        "mag_limit": mag_limit,
        "ref_filter": REF_FILTER,
    })
    return tbl


@add_function_call_str
def omega_cen(fov=60.0,
              imbh_mass=IMBH_MASS,
              time=None,
              epoch=DEFAULT_EPOCH,
              mag_limit=None,
              mass=CLUSTER_MASS,
              distance=DISTANCE,
              r_core=R_CORE,
              r_tide=R_TIDE,
              cusp_slope=0.0,
              turnoff_mass=TURNOFF_MASS,
              turnoff_absmag=TURNOFF_ABSMAG,
              min_mass=MIN_MASS,
              giants=True,
              mass_to_light=MASS_TO_LIGHT,
              metal_poor=True,
              max_stars=2_000_000,
              seed=None):
    """
    Source of an omega-Centauri-like cluster around a central IMBH.

    Only stars whose *projected* position falls inside the requested field
    are generated, but they are drawn from the full 3D cluster, so the
    line-of-sight depth (and hence the velocity distribution) is that of
    the whole cluster. See the module docstring for the model, its
    validation against the literature and its caveats.

    Parameters
    ----------
    fov : float
        [arcsec] Side length of the square field of view. 60 is 1 arcmin;
        MAVIS images 30, MICADO 50 in its wide mode.
    imbh_mass : float
        [Msun] Mass of the central IMBH. Sets the central velocity
        dispersion and the orbits used for time evolution. Default 1e5,
        well above the Haeberle et al. 2024 lower limit of 8.2e3.
    time : astropy.time.Time, str, datetime, or float MJD, optional
        If given, the cluster is evolved from ``epoch`` to ``time`` with
        :func:`evolve_omega_cen` before being returned. With a fixed
        integer ``seed``, different ``time`` values give successive epochs
        of the same realisation.
    epoch : astropy.time.Time, str, datetime, or float MJD
        Reference epoch at which the sampled positions apply.
        Default ``"J2030.0"``.
    mag_limit : float, optional
        [mag] Drop stars fainter than this apparent V magnitude. The
        default keeps everything down to the brown dwarfs; cutting at e.g.
        28 makes ScopeSim much faster.
    mass : float
        [Msun] Total cluster mass (Baumgardt & Hilker 2018: 3.55e6).
    distance : float
        [pc] Distance (Baumgardt & Vasiliev 2021: 5430).
    r_core, r_tide : float
        [pc] Core and tidal radius of the density profile (Harris: 2.37
        and 57.0 *arcmin*, i.e. 3.74 and 90 pc).
    cusp_slope : float
        Slope of a central density cusp inside the IMBH's radius of
        influence. 0 (default) is the observed flat core; 1.75 is a
        Bahcall & Wolf cusp.
    turnoff_mass, turnoff_absmag : float
        [Msun], [mag] Main-sequence turn-off mass and its observed M_V.
        Stars above ``turnoff_mass`` are dropped as dark remnants.
    min_mass : float
        [Msun] Lower end of the IMF.
    giants : bool
        Add the evolved (RGB) component, which dominates the light.
    mass_to_light : float
        [Msun/Lsun] V-band mass-to-light ratio; fixes the number of giants.
    metal_poor : bool
        Prefer the weak-lined Pickles templates ([Fe/H] ~ -1.6).
    max_stars : int
        Guard against generating an unusably large table.
    seed : int, np.random.Generator, or None
        Random seed. ``None`` (default) gives a different cluster on every
        call. At a fixed integer seed, changing ``imbh_mass`` changes only
        the velocities -- positions and magnitudes stay identical.

    Returns
    -------
    scopesim.Source
        One table field. Columns: ``x``, ``y`` [arcsec offsets from the
        cluster centre], ``ref``, ``weight`` (templates are 0 mag in V,
        so ``weight = 10**(-0.4 V)``), ``mag``/``absmag`` [V mag], ``mass``,
        ``phot_mass`` [Msun], ``phase`` (``"MS"``/``"RGB"``),
        ``spec_types``, the 3D state ``x_pc``, ``y_pc``, ``z_pc`` [pc] and
        ``vx``, ``vy``, ``vz`` [km/s], and ``pm_x``, ``pm_y`` [mas/yr].
        Model parameters are in ``src.meta``.

    See Also
    --------
    evolve_omega_cen : move the same stars to another epoch.
    globular_cluster : a simpler, generic cluster with an IMBH.

    Examples
    --------
    >>> from scopesim_templates.stellar import omega_cen, evolve_omega_cen
    >>> src = omega_cen(fov=10, imbh_mass=1e5, mag_limit=26, seed=42)
    >>> later = evolve_omega_cen(src, "J2045")
    """
    tbl = _omega_cen_catalogue(
        fov=fov, imbh_mass=imbh_mass, mass=mass, distance=distance,
        r_core=r_core, r_tide=r_tide, cusp_slope=cusp_slope,
        turnoff_mass=turnoff_mass, turnoff_absmag=turnoff_absmag,
        min_mass=min_mass, giants=giants, mass_to_light=mass_to_light,
        metal_poor=metal_poor, mag_limit=mag_limit, max_stars=max_stars,
        seed=seed)
    if not len(tbl):
        raise ValueError("No stars left in the field; relax mag_limit.")

    # Template per star, then one spectrum per template at 0 mag in V
    pickles_lib = pyckles.SpectralLibrary("pickles", return_style="synphot")
    spec_types = np.asarray(tbl["spec_types"], dtype=object)
    is_ms = tbl["phase"] == "MS"
    if is_ms.any():
        spec_types[is_ms] = _assign_dwarf_templates(
            np.asarray(tbl["phot_mass"][is_ms]),
            list(pickles_lib.table["name"]), metal_poor)
    unique = sorted(set(spec_types))
    spectra = [Spextrum(modelclass=pickles_lib[name]).scale_to_magnitude(
                   filter_curve=REF_FILTER, amplitude=0 * u.mag)
               for name in unique]
    lookup = {name: i for i, name in enumerate(unique)}

    tbl["spec_types"] = np.asarray(spec_types, dtype=str)
    tbl["ref"] = np.array([lookup[name] for name in spec_types], dtype=int)
    tbl["weight"] = 10 ** (-0.4 * tbl["mag"].value)

    meta = dict(tbl.meta)
    tbl.meta.clear()
    src = Source(spectra=spectra, table=tbl)
    src.meta.update(meta)
    src.meta.update({
        "description": "omega Cen (toy model) with central IMBH",
        "ra": RA_OMEGA_CEN,
        "dec": DEC_OMEGA_CEN,
        "epoch": str(epoch),
        "epoch_jyear": to_jyear(epoch),
        "seed": seed_for_meta(seed),
    })
    if time is not None:
        src = evolve_omega_cen(src, time)
    return src


def evolve_omega_cen(source, time, imbh_mass=None):
    """
    Move every star of an :func:`omega_cen` Source to a new epoch.

    Each star is propagated on a Keplerian orbit in the point-mass
    potential of the IMBH *plus* the cluster mass enclosed at its initial
    radius. That is exact well inside the IMBH's radius of influence (about
    90 arcsec for 1e5 Msun, i.e. any realistic field) and, further out, a
    good approximation for epoch differences short compared with the
    crossing time (~1 Myr). Being an orbit solver rather than a linear
    extrapolation, it captures the curvature of the innermost orbits.

    Parameters
    ----------
    source : scopesim.Source
        As returned by :func:`omega_cen` (or by this function).
    time : astropy.time.Time, str, datetime, or float MJD
        Epoch to evolve to (e.g. ``"J2045"``); measured on the same scale
        as ``source.meta["epoch_jyear"]``.
    imbh_mass : float, optional
        [Msun] Black hole mass for the propagation; defaults to the one the
        cluster was built with. Over decades this barely moves anything
        (displacements are ``v dt`` to ~1 per cent); the observable
        signature of the IMBH is in the *velocities*, which are set when
        the cluster is built. To compare masses, rebuild with
        ``omega_cen(imbh_mass=..., seed=<same int>)``.

    Returns
    -------
    scopesim.Source
        A new Source with the same stars and spectra and updated positions,
        velocities, proper motions and epoch. The input is not modified.

    Examples
    --------
    >>> src = omega_cen(fov=10, mag_limit=26, seed=1)
    >>> later = evolve_omega_cen(src, "J2050")
    >>> shift = later.fields[0].field["x"] - src.fields[0].field["x"]
    """
    if len(source.fields) != 1:
        raise ValueError("evolve_omega_cen needs a single-field Source as "
                         "returned by omega_cen().")
    field = source.fields[0]
    meta = dict(field.meta)
    missing = [c for c in _STATE_COLUMNS if c not in field.field.colnames]
    if missing or "epoch_jyear" not in meta:
        raise ValueError(f"Source lacks the omega_cen state ({missing}); "
                         "was it built by omega_cen()?")

    tbl = field.field.copy()
    if imbh_mass is None:
        imbh_mass = meta["imbh_mass"]
    jyear = to_jyear(time)
    dt = jyear - meta["epoch_jyear"]

    if abs(dt) > 1e5:
        warnings.warn(
            f"Evolving by {dt:.3g} yr, a good fraction of the cluster "
            "crossing time. Stars outside the IMBH's radius of influence do "
            "not really follow Keplerian orbits about the centre that long.",
            UserWarning, stacklevel=2)

    pos = np.stack([tbl[c].quantity.to_value(u.pc)
                    for c in ("x_pc", "y_pc", "z_pc")], axis=1)
    vel = np.stack([tbl[c].quantity.to_value(u.km / u.s)
                    for c in ("vx", "vy", "vz")], axis=1) / KMS_PER_PCYR

    if dt:
        potential = _ClusterPotential(meta["cluster_mass"], meta["r_core_pc"],
                                      meta["r_tide_pc"], imbh_mass,
                                      meta["cusp_slope"])
        mu = G_PCYR * potential.central_mass(np.linalg.norm(pos, axis=1))
        pos, vel = kepler_propagate(pos, vel, mu, dt)

    distance = meta["distance_pc"]
    to_mas_yr = KMS_PER_PCYR / (4.740470446 * distance / 1000.0)
    tbl["x"] = pos[:, 0] / distance * RAD2ARCSEC * u.arcsec
    tbl["y"] = pos[:, 1] / distance * RAD2ARCSEC * u.arcsec
    for col, values in zip(("x_pc", "y_pc", "z_pc"), pos.T):
        tbl[col] = values * u.pc
    for col, values in zip(("vx", "vy", "vz"), vel.T):
        tbl[col] = values * KMS_PER_PCYR * u.km / u.s
    tbl["pm_x"] = vel[:, 0] * to_mas_yr * u.mas / u.yr
    tbl["pm_y"] = vel[:, 1] * to_mas_yr * u.mas / u.yr

    spectra = field.spectra
    keys = sorted(spectra)
    if keys != list(range(len(keys))):
        # make refs contiguous so they index the spectra list
        remap = {k: i for i, k in enumerate(keys)}
        tbl["ref"] = [remap[int(r)] for r in tbl["ref"]]
    new = Source(spectra=[spectra[k] for k in keys], table=tbl)
    new.meta.update(meta)
    new.meta.update({
        "epoch": str(time),
        "epoch_jyear": jyear,
        "imbh_mass": float(imbh_mass),
        "evolved_from_jyear": meta["epoch_jyear"],
    })
    return new
