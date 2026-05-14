# -*- coding: utf-8 -*-
"""Globular-cluster star-field template with a central intermediate-mass black hole."""

import numpy as np
from astropy import units as u
from astropy.table import Table

import pyckles
from spextra import Spextrum

from ..rc import Source
from ..utils.general_utils import add_function_call_str
from . import _gillessen
from . import cluster_utils as cu
from . import imf as _imf


__all__ = ["globular_cluster"]


def _sample_salpeter_masses(n_stars, mass_limits, seed):
    """Sample exactly ``n_stars`` masses from a Salpeter (alpha=2.3) IMF."""
    salpeter = _imf.IMF_broken_powerlaw(
        mass_limits=np.array(mass_limits, dtype=float),
        powers=np.array([-2.3]),
    )
    rng = np.random.default_rng(seed)
    mean_m = float(np.mean(mass_limits))  # rough; just to seed total_mass
    total = max(n_stars * mean_m, mass_limits[1] * 2.0)
    # generate_cluster uses np.random; reseed legacy global RNG for determinism
    if seed is not None:
        np.random.seed(int(seed))
    masses = np.array([], dtype=float)
    while masses.size < n_stars:
        new, _, _, _ = salpeter.generate_cluster(total, seed=None)
        masses = np.concatenate([masses, new])
        total *= 1.5
    # Permute then trim to exactly n_stars, so the kept sample is unbiased.
    masses = masses[rng.permutation(masses.size)][:n_stars]
    return masses


def _draw_orbits(n_stars, imbh_mass, a_min_au, a_max_au, eccentricity_cap,
                 rng):
    """Sample random closed orbits and return per-star line-of-sight RV in km/s."""
    log_a = rng.uniform(np.log(a_min_au), np.log(a_max_au), n_stars)
    a_au = np.exp(log_a)

    # Thermal eccentricity, capped to keep orbits well-bound and closed.
    e = np.sqrt(rng.uniform(0.0, 1.0, n_stars))
    e = np.clip(e, 0.0, eccentricity_cap)

    inc = np.arccos(rng.uniform(-1.0, 1.0, n_stars))
    Omega = rng.uniform(0.0, 2.0 * np.pi, n_stars)
    omega = rng.uniform(0.0, 2.0 * np.pi, n_stars)
    M_anom = rng.uniform(0.0, 2.0 * np.pi, n_stars)

    # Mean motion in rad / yr (Kepler 3rd law in solar/AU/yr units:
    # 4 pi^2 a^3 = G M T^2  =>  n^2 = G M / a^3  =>  n = 2 pi sqrt(M / a^3))
    n_rad_yr = 2.0 * np.pi * np.sqrt(imbh_mass / a_au ** 3)

    E = _gillessen._solve_kepler(M_anom, e)
    _, _, vx_orb, vy_orb = _gillessen._orbital_state(a_au, e, n_rad_yr, E)
    _, _, _, _, C, H = _gillessen._thiele_innes(Omega, inc, omega)
    vz_au_yr = C * vx_orb + H * vy_orb
    rv = (vz_au_yr * u.AU / u.yr).to(u.km / u.s).value
    return rv


@add_function_call_str
def globular_cluster(density,
                     fov,
                     distance_modulus,
                     imbh_mass=1e4,
                     filter_name="Generic/Johnson.V",
                     mass_limits=(0.1, 8.0),
                     a_range=(100.0, None),
                     eccentricity_cap=0.95,
                     seed=None):
    """
    Source of a fictional globular cluster with a central IMBH.

    Stars are drawn from a Salpeter IMF, placed uniformly at random in
    a square field of view, and assigned a closed orbit around the
    IMBH at the centre. The line-of-sight velocity of each star at a
    random orbital phase is used to RV-shift its spectrum. Omega
    Centauri is the rough prototype.

    Parameters
    ----------
    density : float
        Stellar surface density [stars / arcsec^2].
    fov : float
        Side length of the square field of view [arcsec]. The total
        number of stars is ``round(density * fov ** 2)``.
    distance_modulus : float
        Apparent magnitudes are computed as ``m = M + distance_modulus``.
        Also implicitly sets the cluster distance
        (``d_pc = 10 ** ((distance_modulus + 5) / 5)``), which scales
        the default outer semi-major axis to the projected field of view.
    imbh_mass : float
        Mass of the central intermediate-mass black hole [solar masses].
        Default 1e4 (Omega Cen-like).
    filter_name : str
        SVO filter ID used to scale the base spectra. Default
        ``"Generic/Johnson.V"`` (matches the mamajek absolute V magnitudes).
    mass_limits : tuple of float
        ``(m_min, m_max)`` in solar masses for the Salpeter sampler.
    a_range : tuple
        ``(a_min, a_max)`` in AU for the log-uniform semi-major axis
        distribution. If ``a_max`` is ``None``, it defaults to half the
        projected field of view at the cluster distance, in AU.
    eccentricity_cap : float
        Upper bound on the (thermal-distribution) eccentricities, kept
        below 1 so all orbits are closed.
    seed : int, optional
        Reproducibility seed.

    Returns
    -------
    scopesim.Source

    Examples
    --------
    >>> from scopesim_templates.stellar import globular_cluster
    >>> src = globular_cluster(density=0.5, fov=20.0,
    ...                        distance_modulus=13.6, seed=42)
    """
    # 1. N stars from density and FOV
    n_stars = int(round(float(density) * float(fov) ** 2))
    if n_stars < 1:
        raise ValueError(
            f"globular_cluster: density * fov**2 = {density * fov**2} "
            "yields zero stars."
        )

    rng = np.random.default_rng(seed)

    # 2. Salpeter masses
    masses = _sample_salpeter_masses(n_stars, mass_limits, seed)

    # 3. spec_types and absolute V mag from mamajek
    spec_types = list(cu.mass2spt(masses))
    Mv = np.asarray(cu.mass2absmag(masses), dtype=float)

    # 4. Base spectra at 0 V mag (one per unique spectral type)
    unique_spts = sorted(set(spec_types))
    pickles = pyckles.SpectralLibrary("pickles", return_style="synphot")
    base = {
        spt: Spextrum(modelclass=pickles[spt]).scale_to_magnitude(
            filter_curve=filter_name, amplitude=0 * u.mag)
        for spt in unique_spts
    }

    # 5. Random closed orbits around the IMBH; line-of-sight velocity
    distance_pc = 10.0 ** ((float(distance_modulus) + 5.0) / 5.0)
    a_min_au = float(a_range[0])
    a_max_au = (float(a_range[1])
                if a_range[1] is not None
                else 0.5 * float(fov) * distance_pc)  # arcsec * pc == AU
    if a_max_au <= a_min_au:
        raise ValueError(
            f"globular_cluster: a_max ({a_max_au}) must exceed a_min "
            f"({a_min_au}). Increase fov, distance_modulus, or a_range[1]."
        )

    rv_kms = _draw_orbits(n_stars, imbh_mass, a_min_au, a_max_au,
                          eccentricity_cap, rng)

    # 6. Per-star unique RV-shifted spectra
    spectra = [base[spec_types[i]].redshift(vel=float(rv_kms[i]))
               for i in range(n_stars)]

    # 7. Apparent magnitudes -> weights
    apparent_V = Mv + float(distance_modulus)
    weight = 10.0 ** (-0.4 * apparent_V)

    # 8. Uniform positions inside the square FOV (orbits decoupled from x/y)
    x = rng.uniform(-fov / 2.0, fov / 2.0, n_stars) * u.arcsec
    y = rng.uniform(-fov / 2.0, fov / 2.0, n_stars) * u.arcsec

    # 9. Assemble Source
    ref = np.arange(n_stars, dtype=int)
    tbl = Table(
        names=["x", "y", "ref", "weight", "mass", "spec_types"],
        data=[x, y, ref, weight, masses, spec_types],
        units=[u.arcsec, u.arcsec, None, None, u.solMass, None],
    )
    src = Source(spectra=spectra, table=tbl)
    src.meta.update({
        "object": "Globular cluster + IMBH (Omega Cen-like)",
        "imbh_mass": imbh_mass,
        "distance_modulus": distance_modulus,
        "distance_pc": distance_pc,
        "filter_name": filter_name,
    })
    return src
