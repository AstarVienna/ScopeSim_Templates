# -*- coding: utf-8 -*-
"""Flux-normalised Sersic image rendering (private helpers).

Ported from the MAVIS 2026 workshop ``hudf.py`` (``sersic_image`` and
friends). The profile is evaluated with numpy directly rather than through
``astropy.modeling.models.Sersic2D`` to avoid the per-call model overhead,
but follows the same parametrisation: ``r_eff`` is the half-light radius
along the major axis and the elliptical radius is
``sqrt(x_maj**2 + (x_min / q)**2)``.

Angles here are *array-frame* angles: degrees counter-clockwise from the
+x (column) axis towards the +y (row) axis. Converting a sky position angle
into this frame is the caller's job.
"""

import numpy as np
from scipy.special import gamma, gammainc, gammaincinv


def sersic_bn(n):
    """Sersic ``b_n``, defined by ``gamma(2n, b_n) = Gamma(2n) / 2``."""
    return gammaincinv(2.0 * np.asarray(n, dtype=float), 0.5)


def sersic_enclosed_fraction(n, n_reff):
    """Fraction of a Sersic profile's light within ``n_reff * r_eff``."""
    return gammainc(2.0 * n, sersic_bn(n) * n_reff ** (1.0 / n))


def sersic_total_flux(r_eff, n, q, amplitude=1.0):
    """Analytic total flux of an elliptical Sersic profile.

    ``amplitude`` is the surface brightness at ``r_eff`` along the major
    axis, matching ``astropy.modeling.models.Sersic2D``. Units are those of
    ``amplitude * r_eff**2``.
    """
    b = sersic_bn(n)
    return (amplitude * 2.0 * np.pi * n * gamma(2.0 * n)
            * np.exp(b) * b ** (-2.0 * n) * r_eff ** 2 * q)


def _sersic_profile(x, y, r_eff, n, q, angle_rad, b_n):
    """Unit-amplitude Sersic surface brightness at (x, y) (same units)."""
    cos_t, sin_t = np.cos(angle_rad), np.sin(angle_rad)
    x_maj = x * cos_t + y * sin_t
    x_min = -x * sin_t + y * cos_t
    radius = np.sqrt(x_maj ** 2 + (x_min / q) ** 2) / r_eff
    return np.exp(-b_n * (radius ** (1.0 / n) - 1.0))


def sersic_image(r_eff, n, q, angle, pixel_scale, n_reff=8.0,
                 oversample=8, core_pixels=24, min_half_size=3,
                 normalise="sum"):
    """Render an elliptical Sersic profile on a square, odd-sized grid.

    The profile is centred on the central pixel and extends to a square of
    half-width ``n_reff * r_eff``. A Sersic profile with large ``n`` is
    cuspy on the scale of a pixel, so the central ``core_pixels`` block is
    integrated over ``oversample**2`` sub-pixel positions instead of being
    sampled once per pixel centre (which overestimates an n=4 core by
    several per cent of the *total* flux).

    Parameters
    ----------
    r_eff : float
        Half-light radius along the major axis, in the same unit as
        `pixel_scale` (normally arcsec).
    n : float
        Sersic index.
    q : float
        Axis ratio b/a, ``0 < q <= 1``.
    angle : float
        Major-axis angle in degrees, counter-clockwise from the array +x
        (column) axis towards +y (row).
    pixel_scale : float
        Size of one pixel, same unit as `r_eff`.
    n_reff : float, optional
        Half-width of the rendered box in units of `r_eff`. Default 8.
    oversample : int, optional
        Sub-sampling factor for the central block. 1 disables it.
    core_pixels : int, optional
        Side length of the over-sampled central block, in pixels.
    min_half_size : int, optional
        Minimum half-width of the image in pixels.
    normalise : {"sum", "total"}, optional
        ``"sum"`` (default) scales the image to unit sum, so that the total
        flux is conserved exactly and the light beyond the box (a few per
        cent for n >= 4) is redistributed inside it. ``"total"`` divides by
        the *analytic* total flux of the untruncated profile, so the surface
        brightness is exact but the image sums to slightly less than one.

    Returns
    -------
    image : np.ndarray
        2D float64 array of odd side length.
    rendered_fraction : float
        Fraction of the analytic total flux contained in the rendered box
        *before* any ``"sum"`` renormalisation (useful bookkeeping).
    """
    if not (r_eff > 0 and n > 0 and 0 < q <= 1 and pixel_scale > 0):
        raise ValueError(f"Invalid Sersic parameters: r_eff={r_eff}, n={n}, "
                         f"q={q}, pixel_scale={pixel_scale}")
    if normalise not in ("sum", "total"):
        raise ValueError(f"normalise must be 'sum' or 'total', "
                         f"got {normalise!r}")

    b_n = float(sersic_bn(n))
    angle_rad = np.deg2rad(angle)
    half = max(int(np.ceil(n_reff * r_eff / pixel_scale)), int(min_half_size))
    size = 2 * half + 1

    axis = (np.arange(size) - half) * pixel_scale
    image = _sersic_profile(axis[None, :], axis[:, None],
                            r_eff, n, q, angle_rad, b_n)

    if oversample > 1 and core_pixels > 0:
        core = min(int(core_pixels), size)
        lo = half - core // 2
        hi = lo + core
        offsets = ((np.arange(oversample) + 0.5) / oversample - 0.5)
        # sub-pixel grid: (core * oversample) samples along each axis
        sub_axis = (axis[lo:hi, None] + offsets[None, :] * pixel_scale)
        sub_axis = sub_axis.ravel()
        sub = _sersic_profile(sub_axis[None, :], sub_axis[:, None],
                              r_eff, n, q, angle_rad, b_n)
        sub = sub.reshape(core, oversample, core, oversample).mean(axis=(1, 3))
        image[lo:hi, lo:hi] = sub

    image *= pixel_scale ** 2 / sersic_total_flux(r_eff, n, q)
    rendered_fraction = float(image.sum())
    if normalise == "sum":
        image /= rendered_fraction
    return image, rendered_fraction
