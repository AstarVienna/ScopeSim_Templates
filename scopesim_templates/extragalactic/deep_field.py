# -*- coding: utf-8 -*-
"""
Hubble Ultra Deep Field analogue, built from a vendored HUDF catalogue.

Every galaxy in the catalogue is rendered as a flux-normalised elliptical
Sersic profile using its published F160W structural parameters, and can
optionally be "clumpified": a mass- and redshift-dependent fraction of the
galaxies have 10-25 per cent of their light moved into compact, power-law
distributed clumps, conserving the total flux. The implementation is a port
of ``hudf.py`` from the MAVIS 2026 workshop feasibility study.

Catalogue
---------
The vendored ``data/hudf_catalogue.fits`` (3152 galaxies within 2.5 arcmin of
the HUDF centre, see ``data/README.md``) combines

* **van der Wel et al. 2012, ApJS 203, 24** (VizieR ``J/ApJS/203/24``):
  GALFIT single-Sersic F160W fits (magnitude, half-light radius, Sersic
  index, axis ratio, position angle), quality flag ``Q = 0``;
* **Skelton et al. 2014, ApJS 214, 24** (VizieR ``J/ApJS/214/24``):
  3D-HST redshifts, FAST stellar masses and F606W photometry.

Please cite both papers when using this template. The catalogue is rebuilt
with ``scripts/make_hudf_catalogue.py`` in the repository (not installed).

Conventions
-----------
* **Positions** are offsets from the field centre (``ra_center``,
  ``dec_center``) in arcsec, ``x = dRA * cos(Dec)`` (positive towards
  **East**) and ``y = dDec`` (positive towards North). They are written to
  each image HDU as ``CRVAL1/2 = x/3600, y/3600`` [deg], which ScopeSim
  interprets as offsets from the optical axis. The absolute centre is kept in
  ``Source.meta`` (``ra_center``, ``dec_center``). Following the ScopeSim
  convention of "x increasing to the right, no RA flip", a detector image
  of this source therefore shows East to the *right*, i.e. it is the mirror
  image of a conventional North-up/East-left sky image.
* **Position angles** in the catalogue are degrees East of North (0 = North,
  90 = East; van der Wel et al. 2012). In the (x, y) frame above, a major
  axis at PA points along ``(sin PA, cos PA)``, i.e. at
  ``theta = 90 deg - PA`` counter-clockwise from +x. This keeps each
  galaxy's orientation consistent with the positions (the whole scene is
  mirrored consistently, not galaxy by galaxy). Note that the workshop
  ``hudf.py`` used ``theta = PA + 90``, which is the convention for an
  East-left frame and mirrors every galaxy relative to the positions.
* **Flux**: each image sums to ``10**(-0.4 * m_F160W)`` (times the fraction
  left after clumps) and its spectrum is normalised to 0 mag(AB) in
  `filter_name`, so ``sum(image) * spectrum`` reproduces the catalogue
  magnitude. With ``conserve_flux=True`` (default) the Sersic wings beyond
  ``n_reff * r_eff`` are folded back into the rendered box.

Clump model (``clumps=True``)
-----------------------------
* Which galaxies are clumpy: Guo et al. 2015, ApJ 800, 39 clumpy fractions vs
  stellar mass and redshift (:data:`CLUMPY_FRACTION`), plus a Sersic-index
  cut ``clump_sersic_max`` (default 2.5) to keep bulges smooth.
* Total clump light fraction: uniform in 10-25 per cent (Foerster Schreiber
  et al. 2011, ApJ 739, 45).
* Clump light function ``dN/df ~ f**-slope``, slope 2 (Dessauges-Zavadsky &
  Adamo 2018, MNRAS 479, L118), individual fractions 0.5-15 per cent.
* Positions trace the galaxy's own Sersic light profile.
* Clumps are point sources with a blue compact starburst SED (default
  ``brown/UM461``), and their flux is taken out of the smooth component so
  the galaxy total is unchanged.

The clump populations are random draws, not measurements of the real HUDF
galaxies.
"""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from spextra import Passband, Spextrum
from synphot import ConstFlux1D, Empirical1D, SourceSpectrum, SpectralElement
from synphot.units import PHOTLAM

from ..rc import Source
from ..utils.general_utils import add_function_call_str
from ._sersic import sersic_bn, sersic_image


__all__ = ["hubble_ultra_deep_field",
           "load_hudf_catalogue",
           "select_field",
           "richest_field",
           "clumpy_fraction",
           "clump_catalogue"]


_DEFAULT_CATALOGUE = Path(__file__).parent / "data" / "hudf_catalogue.fits"

#: HUDF ACS field centre (Beckwith et al. 2006, AJ 132, 1729), degrees.
HUDF_RA = 53.1625
HUDF_DEC = -27.79142

#: Total fraction of galaxy light in clumps, drawn uniformly per galaxy
#: (Foerster Schreiber et al. 2011).
CLUMP_FLUX_FRACTION = (0.10, 0.25)
#: ``dN/df ~ f**-slope`` (Dessauges-Zavadsky & Adamo 2018: inferred 2.0).
CLUMP_SLOPE = 2.0
#: Range of individual clump light fractions (Foerster Schreiber et al. 2011).
CLUMP_FRACTION_RANGE = (0.005, 0.15)
#: Clumpy fraction of galaxies, digitised from Guo et al. 2015. Keys are the
#: upper edges of log10(M*/Msun) bins; values are (z, fraction) knots,
#: linearly interpolated in z and held flat outside.
CLUMPY_FRACTION = {
    9.8:  [(0.5, 0.60), (1.5, 0.60), (3.0, 0.60)],
    10.6: [(0.5, 0.40), (1.5, 0.55), (3.0, 0.55)],
    99.0: [(0.5, 0.15), (1.0, 0.30), (2.0, 0.55), (3.0, 0.55)],
}

# Column names used internally (those of the vendored catalogue), and the
# aliases accepted from user catalogues (e.g. the workshop ``hudf.py`` cache).
_COLUMN_ALIASES = {
    "m_F160W": ("mag",),
    "r_eff_arcsec": ("re", "r_eff"),
    "sersic_n": ("n",),
    "axis_ratio": ("q",),
    "position_angle_deg": ("pa",),
    "log_mstar": ("logm",),
    "m_F606W": ("mag_f606w",),
}
_REQUIRED = ("ra", "dec", "z", "m_F160W", "r_eff_arcsec", "sersic_n",
             "axis_ratio", "position_angle_deg")


# ---------------------------------------------------------------------------
# catalogue handling
# ---------------------------------------------------------------------------

def load_hudf_catalogue(catalogue_path=None, ra_center=HUDF_RA,
                        dec_center=HUDF_DEC):
    """Read an HUDF galaxy catalogue and add field-centre offsets.

    Parameters
    ----------
    catalogue_path : str or Path, optional
        FITS (or any astropy-readable) table. Default: the vendored
        catalogue. Required columns: ``ra, dec, z, m_F160W, r_eff_arcsec,
        sersic_n, axis_ratio, position_angle_deg`` (the workshop names
        ``mag, re, n, q, pa`` are accepted too). Optional: ``id``,
        ``log_mstar`` (``logm``) and ``m_F606W``.
    ra_center, dec_center : float
        Field centre in degrees, from which the offsets are computed.

    Returns
    -------
    astropy.table.Table
        A new table with the canonical column names plus ``x``, ``y``:
        offsets from the centre in arcsec, ``x`` positive towards East,
        ``y`` towards North.
    """
    path = (Path(catalogue_path) if catalogue_path is not None
            else _DEFAULT_CATALOGUE)
    tbl = Table.read(path)          # a fresh object, safe to modify
    for name, aliases in _COLUMN_ALIASES.items():
        if name not in tbl.colnames:
            for alias in aliases:
                if alias in tbl.colnames:
                    tbl.rename_column(alias, name)
                    break

    missing = [c for c in _REQUIRED if c not in tbl.colnames]
    if missing:
        raise KeyError(f"Catalogue {path} lacks required columns {missing}.")
    if "id" not in tbl.colnames:
        tbl["id"] = np.arange(len(tbl))
    for optional in ("log_mstar", "m_F606W"):
        if optional not in tbl.colnames:
            tbl[optional] = np.full(len(tbl), np.nan)

    centre = SkyCoord(ra_center * u.deg, dec_center * u.deg)
    rows = SkyCoord(np.asarray(tbl["ra"], dtype=float) * u.deg,
                    np.asarray(tbl["dec"], dtype=float) * u.deg)
    dra, ddec = centre.spherical_offsets_to(rows)
    tbl["x"] = dra.to_value(u.arcsec)
    tbl["y"] = ddec.to_value(u.arcsec)
    tbl.meta["catalogue"] = str(path)
    tbl.meta["ra_center"] = float(ra_center)
    tbl.meta["dec_center"] = float(dec_center)
    return tbl


def select_field(catalogue, size, x0=0.0, y0=0.0, margin=0.0,
                 mag_limit=None):
    """Cut a square patch out of a catalogue with ``x``/``y`` columns.

    Parameters
    ----------
    catalogue : Table
        Output of :func:`load_hudf_catalogue`.
    size : float or None
        Side length in arcsec. ``None`` keeps the whole catalogue.
    x0, y0 : float, optional
        Patch centre in arcsec, in the catalogue's ``x``/``y`` frame.
    margin : float, optional
        Extra border in arcsec, so that galaxies just outside the patch can
        still spill light into it. Default 0.
    mag_limit : float, optional
        Drop galaxies fainter than this ``m_F160W``.

    Returns
    -------
    Table
        A copy of the selected rows with ``x``/``y`` re-referenced to the
        patch centre.
    """
    keep = np.ones(len(catalogue), dtype=bool)
    if mag_limit is not None:
        keep &= np.asarray(catalogue["m_F160W"]) <= float(mag_limit)
    if size is not None:
        half = 0.5 * float(size) + float(margin)
        keep &= ((np.abs(np.asarray(catalogue["x"]) - x0) <= half)
                 & (np.abs(np.asarray(catalogue["y"]) - y0) <= half))
    out = catalogue[keep].copy()
    out["x"] = np.asarray(out["x"], dtype=float) - x0
    out["y"] = np.asarray(out["y"], dtype=float) - y0
    out.meta.update({"field_size": size, "field_x0": x0, "field_y0": y0,
                     "field_margin": margin})
    return out


def richest_field(catalogue, size, step=2.0, mag_limit=25.0, z_range=None,
                  weight_by_flux=True, max_offset=None):
    """Find the square patch with the most (or brightest) galaxies.

    Scans a grid of trial centres and scores each by the number of galaxies
    passing the cuts, or by their summed F160W flux. The cuts only affect
    the score, not what :func:`select_field` later returns.

    Parameters
    ----------
    catalogue : Table
        Output of :func:`load_hudf_catalogue`.
    size : float
        Side length of the patch in arcsec.
    step : float, optional
        Grid spacing of trial centres in arcsec. Default 2.
    mag_limit : float, optional
        Only galaxies brighter than this count. Default 25.
    z_range : tuple, optional
        ``(z_min, z_max)`` scoring cut, e.g. ``(1, 3)`` to target the epoch
        where clumpy galaxies dominate.
    weight_by_flux : bool, optional
        Score by summed flux (default) instead of galaxy count.
    max_offset : float, optional
        Keep the patch centre within this many arcsec of the catalogue
        centre in each axis.

    Returns
    -------
    tuple
        ``(x0, y0, score)``. Convert to sky coordinates with e.g.
        ``SkyCoord(ra, dec, unit="deg").spherical_offsets_by(x0 * u.arcsec,
        y0 * u.arcsec)`` and pass them as ``ra_center``/``dec_center`` to
        :func:`hubble_ultra_deep_field`.
    """
    x = np.asarray(catalogue["x"], dtype=float)
    y = np.asarray(catalogue["y"], dtype=float)
    mag = np.asarray(catalogue["m_F160W"], dtype=float)
    keep = mag < mag_limit
    if z_range is not None:
        z = np.asarray(catalogue["z"], dtype=float)
        keep &= (z >= z_range[0]) & (z <= z_range[1])
    x, y, mag = x[keep], y[keep], mag[keep]
    weights = 10 ** (-0.4 * mag) if weight_by_flux else np.ones(len(mag))

    reach = max(np.abs(catalogue["x"]).max(),
                np.abs(catalogue["y"]).max()) - size / 2.0
    if max_offset is not None:
        reach = min(reach, max_offset)
    reach = max(reach, 0.0)
    grid = np.arange(-reach, reach + 0.5 * step, step)
    half = size / 2.0

    # (n_grid, n_gal) membership per axis; score = in_x @ (w * in_y).T
    in_x = (np.abs(x[None, :] - grid[:, None]) < half).astype(float)
    in_y = (np.abs(y[None, :] - grid[:, None]) < half).astype(float)
    scores = (in_x * weights[None, :]) @ in_y.T     # [ix, iy]
    ix, iy = np.unravel_index(np.argmax(scores), scores.shape)
    return float(grid[ix]), float(grid[iy]), float(scores[ix, iy])


# ---------------------------------------------------------------------------
# clumps
# ---------------------------------------------------------------------------

def _pa_to_theta(pa_deg):
    """Sky PA (deg E of N) -> angle CCW from +x (deg), for x = East."""
    return 90.0 - np.asarray(pa_deg, dtype=float)


def clumpy_fraction(logmass, z, table=None):
    """Probability that a galaxy is clumpy, after Guo et al. 2015.

    Parameters
    ----------
    logmass, z : array_like
        log10 stellar mass [Msun] and redshift. A missing (NaN) mass is
        treated as low mass, the most common case.
    table : dict, optional
        Overrides :data:`CLUMPY_FRACTION`.

    Returns
    -------
    np.ndarray
        Clumpy fraction per input, held flat outside the tabulated z range.
    """
    table = CLUMPY_FRACTION if table is None else table
    logmass, z = np.broadcast_arrays(np.atleast_1d(np.asarray(logmass, float)),
                                     np.atleast_1d(np.asarray(z, float)))
    out = np.zeros(logmass.shape)
    edges = sorted(table)
    mass = np.where(np.isfinite(logmass), logmass, edges[0] - 1.0)
    lower = -np.inf
    for edge in edges:
        knots = np.asarray(table[edge], dtype=float)
        in_bin = (mass > lower) & (mass <= edge)
        out[in_bin] = np.interp(z[in_bin], knots[:, 0], knots[:, 1])
        lower = edge
    return out


def _power_law_samples(rng, size, slope, low, high):
    """Draw from ``dN/df ~ f**-slope`` on ``[low, high]``."""
    uniform = rng.random(size)
    if np.isclose(slope, 1.0):
        return low * (high / low) ** uniform
    power = 1.0 - slope
    return (low ** power + uniform * (high ** power - low ** power)) \
        ** (1.0 / power)


def _power_law_mean(slope, low, high):
    """Mean of ``dN/df ~ f**-slope`` truncated to ``[low, high]``."""
    if np.isclose(slope, 1.0):
        return (high - low) / np.log(high / low)
    if np.isclose(slope, 2.0):
        return np.log(high / low) / (1.0 / low - 1.0 / high)
    a, b = 1.0 - slope, 2.0 - slope
    return (a / b) * (high ** b - low ** b) / (high ** a - low ** a)


def _sersic_radii(rng, count, r_eff, n, n_reff, radial_index, n_grid=2048):
    """Semi-major-axis radii drawn from the Sersic light distribution."""
    b = float(sersic_bn(n))
    r = np.linspace(0.0, n_reff * r_eff, n_grid)
    r_safe = np.maximum(r, 1e-12 * r_eff)
    weight = (np.exp(-b * ((r_safe / r_eff) ** (1.0 / n) - 1.0))
              * r_safe ** (1.0 + radial_index))
    cdf = np.cumsum(weight)
    cdf /= cdf[-1]
    return np.interp(rng.random(count), cdf, r)


def clump_catalogue(galaxies, flux_fraction=CLUMP_FLUX_FRACTION,
                    slope=CLUMP_SLOPE, fraction_range=CLUMP_FRACTION_RANGE,
                    radial_index=0.0, sersic_max=2.5, n_reff=8.0,
                    clumpy_table=None, mag_limit=None, rescale=True,
                    max_total=0.6, seed=None):
    """Sample a population of point-like clumps for each galaxy.

    See the module docstring for the literature behind every default.

    Parameters
    ----------
    galaxies : Table
        Output of :func:`load_hudf_catalogue` / :func:`select_field`.
    flux_fraction : tuple or float, optional
        Total fraction of a clumpy galaxy's light placed in clumps. A 2-tuple
        is drawn uniformly per galaxy.
    slope : float, optional
        Power-law slope of the clump light-fraction distribution.
    fraction_range : tuple, optional
        ``(min, max)`` light fraction of an individual clump.
    radial_index : float, optional
        Extra ``r**radial_index`` weighting of clump radii (0: trace light).
    sersic_max : float or None, optional
        Only galaxies with ``n <= sersic_max`` may be clumpy.
    n_reff : float, optional
        Clumps are placed within this many effective radii.
    clumpy_table : dict, optional
        Overrides :data:`CLUMPY_FRACTION`.
    mag_limit : float, optional
        Only galaxies brighter than this may be clumpy.
    rescale : bool, optional
        If True (default) rescale each galaxy's clump fractions to sum to
        exactly the drawn total; if False keep the raw power-law draws
        (capped at `max_total`).
    max_total : float, optional
        Upper limit of the summed clump fraction with ``rescale=False``.
    seed : int, np.random.Generator or None, optional
        Random state. ``None`` (default) is not reproducible.

    Returns
    -------
    Table
        One row per clump: ``host`` (row index into `galaxies`), ``host_id``,
        ``x``, ``y`` (arcsec, frame of `galaxies`), ``dx``, ``dy`` (from the
        host centre), ``frac`` (fraction of host flux), ``mag`` (F160W AB)
        and ``z``. ``meta["taken"]`` holds, per galaxy, the summed fraction
        moved into clumps.
    """
    rng = np.random.default_rng(seed)
    low, high = fraction_range
    mean_fraction = _power_law_mean(slope, low, high)

    n_gal = len(galaxies)
    sersic_n = np.asarray(galaxies["sersic_n"], dtype=float)
    eligible = np.ones(n_gal, dtype=bool)
    if sersic_max is not None:
        eligible &= sersic_n <= sersic_max
    if mag_limit is not None:
        eligible &= np.asarray(galaxies["m_F160W"]) <= mag_limit
    probability = clumpy_fraction(galaxies["log_mstar"], galaxies["z"],
                                  table=clumpy_table)
    is_clumpy = eligible & (rng.random(n_gal) < probability)

    cols = {k: [] for k in ("host", "host_id", "x", "y", "dx", "dy", "frac",
                            "mag", "z")}
    for index in np.flatnonzero(is_clumpy):
        gal = galaxies[index]
        total = (rng.uniform(*flux_fraction) if np.iterable(flux_fraction)
                 else float(flux_fraction))
        expected = max(int(round(total / mean_fraction)), 1)
        count = max(int(rng.poisson(expected)), 1)
        fractions = _power_law_samples(rng, count, slope, low, high)
        if rescale:
            fractions *= total / fractions.sum()
        else:
            for _ in range(50):
                if fractions.sum() <= max_total:
                    break
                fractions = _power_law_samples(rng, count, slope, low, high)
            else:
                fractions *= max_total / fractions.sum()

        q = float(gal["axis_ratio"])
        radius = _sersic_radii(rng, count, float(gal["r_eff_arcsec"]),
                               float(gal["sersic_n"]), n_reff, radial_index)
        azimuth = rng.uniform(0.0, 2.0 * np.pi, count)
        major = radius * np.cos(azimuth)
        minor = radius * q * np.sin(azimuth)
        angle = np.deg2rad(_pa_to_theta(gal["position_angle_deg"]))
        dx = major * np.cos(angle) - minor * np.sin(angle)
        dy = major * np.sin(angle) + minor * np.cos(angle)

        cols["host"].extend([index] * count)
        cols["host_id"].extend([int(gal["id"])] * count)
        cols["x"].extend(float(gal["x"]) + dx)
        cols["y"].extend(float(gal["y"]) + dy)
        cols["dx"].extend(dx)
        cols["dy"].extend(dy)
        cols["frac"].extend(fractions)
        cols["mag"].extend(float(gal["m_F160W"]) - 2.5 * np.log10(fractions))
        cols["z"].extend([float(gal["z"])] * count)

    dtypes = {"host": int, "host_id": int}
    clumps = Table({k: np.asarray(v, dtype=dtypes.get(k, float))
                    for k, v in cols.items()})
    taken = np.zeros(n_gal)
    if len(clumps):
        np.add.at(taken, clumps["host"], np.asarray(clumps["frac"]))
    clumps.meta.update({
        "slope": slope, "fraction_range": tuple(fraction_range),
        "radial_index": radial_index, "sersic_max": sersic_max,
        "n_clumpy": int(is_clumpy.sum()), "n_galaxies": n_gal,
        "rescale": bool(rescale), "taken": taken,
    })
    return clumps


# ---------------------------------------------------------------------------
# Source assembly
# ---------------------------------------------------------------------------

class _SpectrumCache:
    """Redshifted SEDs normalised to 0 mag(AB) in a filter, keyed by z-bin.

    ``Spextrum.redshift().scale_to_magnitude()`` costs ~0.2 s per call
    (it re-reads the filter every time), which dominated the build time of
    large fields. Here the filter is read once, the template is redshifted
    by stretching its wavelength grid, and the AB normalisation is done by
    direct photon-counting integration::

        10**(-0.4 m_AB) = int(photlam * T dlam) / int(photlam_AB0 * T dlam)

    which is the synphot definition of an AB magnitude.
    """

    def __init__(self, filter_name, z_step):
        self.filter_name = filter_name
        self.z_step = z_step
        self._base = {}
        self._keys = {}
        self.spectra = {}
        self._filter = None

    def _filter_grid(self):
        if self._filter is None:
            passband = (self.filter_name
                        if isinstance(self.filter_name, SpectralElement)
                        else Passband(self.filter_name))
            wave = passband.waveset.to_value(u.AA)
            trans = np.asarray(passband(wave * u.AA), dtype=float)
            ab0 = SourceSpectrum(ConstFlux1D, amplitude=0 * u.ABmag)
            ab0_counts = np.trapezoid(
                ab0(wave * u.AA).to_value(PHOTLAM) * trans, wave)
            self._filter = (wave, trans, ab0_counts)
        return self._filter

    def _base_spectrum(self, sed):
        name = sed if isinstance(sed, str) else id(sed)
        if name not in self._base:
            spec = Spextrum(sed) if isinstance(sed, str) else sed
            if spec.waveset is not None:
                wave = spec.waveset.to_value(u.AA)
            else:                   # analytic model, e.g. a flat spectrum
                wave = np.geomspace(100.0, 3e5, 4000)
            flux = spec(wave * u.AA).to_value(PHOTLAM)
            self._base[name] = (wave, flux)
        return name, self._base[name]

    def key(self, sed, z):
        """Integer spectrum ref for `sed` at redshift `z` (binned)."""
        z_bin = (float(np.round(z / self.z_step) * self.z_step)
                 if self.z_step else float(z))
        z_bin = max(z_bin, 0.0)
        name, (wave, flux) = self._base_spectrum(sed)
        if (name, z_bin) not in self._keys:
            f_wave, f_trans, ab0_counts = self._filter_grid()
            wave_z = wave * (1.0 + z_bin)
            inside = (wave_z > f_wave[0]) & (wave_z < f_wave[-1])
            grid = np.union1d(f_wave, wave_z[inside])
            counts = np.trapezoid(
                np.interp(grid, wave_z, flux, left=0.0, right=0.0)
                * np.interp(grid, f_wave, f_trans), grid)
            if not counts > 0:
                raise ValueError(
                    f"SED {sed!r} at z={z_bin} has no flux in "
                    f"{self.filter_name}; cannot normalise it.")
            spec = SourceSpectrum(Empirical1D, points=wave_z * u.AA,
                                  lookup_table=flux * (ab0_counts / counts)
                                  * PHOTLAM)
            ref = len(self._keys)
            self._keys[(name, z_bin)] = ref
            self.spectra[ref] = spec
        return self._keys[(name, z_bin)]


def _image_hdu(image, x, y, pixel_scale):
    """ImageHDU whose WCS puts the image centre at offset (x, y) [arcsec]."""
    ny, nx = image.shape
    hdr = fits.Header()
    hdr["CTYPE1"], hdr["CTYPE2"] = "RA---TAN", "DEC--TAN"
    hdr["CUNIT1"], hdr["CUNIT2"] = "deg", "deg"
    hdr["CRPIX1"] = (nx + 1) / 2
    hdr["CRPIX2"] = (ny + 1) / 2
    hdr["CRVAL1"] = x / 3600.0
    hdr["CRVAL2"] = y / 3600.0
    # Negative CDELT1 is the FITS sky convention; ScopeSim ignores its sign
    # when placing image fields, so the array is never mirrored.
    hdr["CDELT1"] = -pixel_scale / 3600.0
    hdr["CDELT2"] = pixel_scale / 3600.0
    hdr["BUNIT"] = "photlam"     # binned weights; flux = sum(img) * spectrum
    return fits.ImageHDU(data=image.astype(np.float32), header=hdr)


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
                            catalogue_path=None,
                            *,
                            margin=0.0,
                            n_reff=8.0,
                            max_pixels=1024,
                            conserve_flux=True,
                            z_step=0.05,
                            clumps=False,
                            clump_flux_fraction=CLUMP_FLUX_FRACTION,
                            clump_slope=CLUMP_SLOPE,
                            clump_fraction_range=CLUMP_FRACTION_RANGE,
                            clump_sersic_max=2.5,
                            clump_radial_index=0.0,
                            clump_sed="brown/UM461",
                            seed=None):
    """
    Source of a Hubble-Ultra-Deep-Field galaxy population.

    Loads a vendored catalogue of HUDF galaxies (van der Wel+2012 F160W
    Sersic fits cross-matched with Skelton+2014 3D-HST redshifts and
    masses), trims it to a square FOV and magnitude limit, and renders each
    galaxy as a flux-normalised elliptical Sersic image. Optionally, a
    fraction of the galaxies is "clumpified" (see ``clumps``). See the
    module docstring for the position, position-angle and flux conventions.

    Each galaxy gets one of two spextra SEDs based on a Sersic-index cut
    (``n > n_cut`` is treated as early-type), redshifted to the catalogue
    ``z`` (binned by `z_step`) and scaled so the galaxy has its catalogue
    ``m_F160W`` in `filter_name`.

    Parameters
    ----------
    pixel_scale : float
        Pixel scale of the rendered galaxy images in arcsec/pixel. Default
        0.06. ScopeSim resamples onto the detector conserving flux, so this
        only needs to be fine enough to resolve the galaxies' structure.
    fov : float, optional
        Side length of the square FOV in arcsec, centred on
        (`ra_center`, `dec_center`). Galaxies whose centres fall inside
        (plus `margin`) are kept. ``None`` (default) keeps the whole
        catalogue (2.5 arcmin radius around the HUDF centre).
    ra_center, dec_center : float
        FOV centre in degrees. Default is the HUDF centre (Beckwith+2006).
        Galaxy positions are offsets from this point (see module docstring).
    mag_limit : float, optional
        Drop galaxies fainter than this F160W AB magnitude. ``None`` keeps
        all (the vendored catalogue reaches ~29 mag).
    filter_name : str or synphot.SpectralElement
        SVO filter ID (or a synphot bandpass) in which the catalogue
        magnitude is applied. Default ``"HST/WFC3_IR.F160W"``, the band of
        the GALFIT magnitudes.
    sed_early, sed_late : str, spextra.Spextrum or synphot.SourceSpectrum
        SEDs for ``n > n_cut`` and ``n <= n_cut`` galaxies. Defaults are
        Brown+2014 archetypes ``NGC0584`` (E0) and ``NGC4254`` (Sc).
    n_cut : float
        Sersic-index threshold between early- and late-type SEDs.
    catalogue_path : str or Path, optional
        Use another catalogue instead of the vendored one. See
        :func:`load_hudf_catalogue` for the required columns.

    Other Parameters
    ----------------
    margin : float
        Extra border in arcsec around `fov` within which galaxies are still
        included (their light may spill into the FOV). Default 0.
    n_reff : float
        Each galaxy is rendered out to ``n_reff * r_eff``. Default 8.
    max_pixels : int or None
        Maximum side length of a galaxy image. Galaxies that would be larger
        are rendered on a coarser grid (not truncated). Default 1024.
    conserve_flux : bool
        If True (default), each rendered profile is renormalised to unit sum
        so the total flux matches the catalogue exactly. If False, the
        profile keeps its exact analytic surface brightness and the light
        beyond ``n_reff * r_eff`` (~3.6 % for n = 4, n_reff = 8) is lost.
    z_step : float
        Redshift bin width for sharing spectra between galaxies. Default
        0.05 (1.7 % error in 1+z at z = 2). 0 disables binning.
    clumps : bool
        If True, move part of the light of clumpy galaxies into point-source
        clumps (see module docstring). Default False.
    clump_flux_fraction : tuple or float
        Total light fraction in clumps per clumpy galaxy. Default
        (0.10, 0.25), drawn uniformly.
    clump_slope : float
        Slope of ``dN/df ~ f**-slope``. Default 2.
    clump_fraction_range : tuple
        Light-fraction range of an individual clump. Default (0.005, 0.15).
    clump_sersic_max : float or None
        Only galaxies with ``n <= clump_sersic_max`` may be clumpy.
    clump_radial_index : float
        Extra ``r**index`` weighting of clump radii (0: clumps trace light).
    clump_sed : str or spextra.Spextrum
        SED of the clumps. Default ``"brown/UM461"`` (blue compact dwarf).
    seed : int, np.random.Generator or None
        Random state for the clump draws. ``None`` (default) is not
        reproducible. The global numpy random state is never touched.

    Returns
    -------
    scopesim.Source
        One image field per galaxy and, with ``clumps=True``, one table
        field holding all clumps.

    Examples
    --------
    >>> from scopesim_templates.extragalactic import hubble_ultra_deep_field
    >>> src = hubble_ultra_deep_field(fov=20.0, mag_limit=27.0)
    >>> src = hubble_ultra_deep_field(fov=10.0, pixel_scale=0.004,
    ...                               clumps=True, seed=42)
    """
    cat = load_hudf_catalogue(catalogue_path, ra_center, dec_center)
    sel = select_field(cat, fov, margin=margin, mag_limit=mag_limit)
    if len(sel) == 0:
        raise ValueError(
            f"hubble_ultra_deep_field: zero galaxies after fov={fov} and "
            f"mag_limit={mag_limit} cuts on {cat.meta['catalogue']}.")

    clump_tbl = None
    taken = np.zeros(len(sel))
    if clumps:
        clump_tbl = clump_catalogue(
            sel, flux_fraction=clump_flux_fraction, slope=clump_slope,
            fraction_range=clump_fraction_range,
            radial_index=clump_radial_index, sersic_max=clump_sersic_max,
            n_reff=n_reff, seed=seed)
        taken = clump_tbl.meta["taken"]

    cache = _SpectrumCache(filter_name, z_step)
    src = Source()
    rendered = np.empty(len(sel))
    n_early = 0
    for i, row in enumerate(sel):
        r_eff = float(row["r_eff_arcsec"])
        n = float(row["sersic_n"])
        scale = float(pixel_scale)
        if max_pixels is not None:
            scale = max(scale, 2.0 * n_reff * r_eff / max_pixels)
        image, rendered[i] = sersic_image(
            r_eff, n, float(row["axis_ratio"]),
            float(_pa_to_theta(row["position_angle_deg"])),
            pixel_scale=scale, n_reff=n_reff,
            normalise="sum" if conserve_flux else "total")
        image *= 10 ** (-0.4 * float(row["m_F160W"])) * (1.0 - taken[i])

        early = n > n_cut
        n_early += early
        ref = cache.key(sed_early if early else sed_late, float(row["z"]))
        hdu = _image_hdu(image, float(row["x"]), float(row["y"]), scale)
        hdu.header["SPEC_REF"] = ref
        hdu.header["GAL_ID"] = int(row["id"])
        # Build each field through the public Source constructor, then
        # collect the fields in one Source. ``Source.__add__`` deep-copies
        # the whole accumulated Source each time, which is O(n^2).
        fld = Source(image_hdu=hdu, spectra=[cache.spectra[ref]]).fields[0]
        fld.field.header["SPEC_REF"] = ref
        fld.spectra = {ref: cache.spectra[ref]}
        fld.meta.update({"galaxy_id": int(row["id"]), "z": float(row["z"]),
                         "m_F160W": float(row["m_F160W"])})
        src.fields.append(fld)

    n_clumps = 0
    if clump_tbl is not None and len(clump_tbl):
        n_clumps = len(clump_tbl)
        refs = [cache.key(clump_sed, z) for z in clump_tbl["z"]]
        spectra = {ref: cache.spectra[ref] for ref in set(refs)}
        clump_src = Source(x=np.asarray(clump_tbl["x"]),
                           y=np.asarray(clump_tbl["y"]),
                           ref=np.zeros(n_clumps, dtype=int),
                           weight=10 ** (-0.4 * np.asarray(clump_tbl["mag"])),
                           spectra=[cache.spectra[refs[0]]])
        fld = clump_src.fields[0]
        fld.field["ref"] = np.asarray(refs, dtype=int)
        fld.field["host_id"] = np.asarray(clump_tbl["host_id"])
        fld.field["frac"] = np.asarray(clump_tbl["frac"])
        fld.spectra = spectra
        src.fields.append(fld)

    src.meta.update({
        "object": "Hubble Ultra Deep Field",
        "ra_center": float(ra_center),
        "dec_center": float(dec_center),
        "fov": fov,
        "margin": margin,
        "mag_limit": mag_limit,
        "pixel_scale": pixel_scale,
        "filter_name": str(filter_name),
        "sed_early": str(sed_early),
        "sed_late": str(sed_late),
        "n_cut": n_cut,
        "n_galaxies": len(sel),
        "n_early": int(n_early),
        "n_late": int(len(sel) - n_early),
        "n_clumps": n_clumps,
        "n_clumpy_galaxies": int(np.count_nonzero(taken)),
        "n_spectra": len(cache.spectra),
        "n_reff": n_reff,
        "conserve_flux": conserve_flux,
        "rendered_fraction_median": float(np.median(rendered)),
        "rendered_fraction_min": float(np.min(rendered)),
        "total_flux_F160W_mag": float(
            -2.5 * np.log10(np.sum(10 ** (-0.4 * np.asarray(sel["m_F160W"]))))),
        "position_convention": "x = dRA*cos(Dec) (East, +right), y = dDec "
                               "[arcsec]; CRVAL = offset from ra/dec_center",
        "pa_convention": "theta = 90 deg - PA (PA east of north)",
        "catalogue": cat.meta["catalogue"],
        "references": "van der Wel+2012 ApJS 203 24; Skelton+2014 ApJS 214 "
                      "24",
        "seed": seed if isinstance(seed, (int, np.integer, type(None)))
                else repr(seed),
    })
    return src
