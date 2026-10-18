# -*- coding: utf-8 -*-
"""
Procedural generator for semi-realistic 2D nebulosity images.

The approach combines a few standard, FFT-based techniques so that a 2k x 2k
image takes about a second on a laptop (4k x 4k well under 10 s):

1. A Gaussian random field with a power-law spectrum, P(k) ~ k^-beta. The
   default beta of 2.7 is close to the Kolmogorov value (8/3) and to the
   slope measured from narrow-band images of NGC 6888.
2. A "ridged" multi-octave component (Musgrave-style ridged multifractal).
   Each octave band g_i is turned into ridges with exp(-p |g_i|), which
   produces thin filaments and sheet edges.
3. Optional domain warping (sampling the field at x + A D(x)), which
   stretches the structure into wispy, wind-blown shapes.
4. A lognormal transform, I = exp(sigma L). This matches the observed
   density PDF of the ISM and gives the strongly skewed brightness
   distribution typical of emission nebulae.
5. An optional large-scale envelope, e.g. a noisy elliptical shell similar
   to the Crescent Nebula.
"""

import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy import fft, ndimage

from ..rc import __config__


__all__ = ["nebulosity"]


def _kgrid(shape, full_shape=None):
    """
    Return |k| on the ``rfft2`` grid for ``shape``.

    Units are cycles per pixel of ``full_shape`` (default: ``shape``), so
    a downsampled grid shares the wavenumbers of the full-resolution one.
    """
    full_shape = full_shape or shape
    ky = fft.fftfreq(shape[0], d=full_shape[0] / shape[0])
    kx = fft.rfftfreq(shape[1], d=full_shape[1] / shape[1])
    return np.hypot(kx[None, :], ky[:, None]).astype(np.float32)


def _normalise(arr):
    """Shift and scale to zero mean and unit standard deviation."""
    arr -= arr.mean()
    std = arr.std()
    if std > 0:
        arr /= std
    return arr


def _irfft2(spec, shape):
    return fft.irfft2(spec, s=shape, workers=-1)


def _white_spectrum(shape, rng):
    noise = rng.standard_normal(shape, dtype=np.float32)
    return fft.rfft2(noise, workers=-1)


def _crop_spectrum(spec, full_shape, shape):
    """Keep only the wavenumbers of ``spec`` that fit on a smaller grid."""
    n_top = (shape[0] + 1) // 2
    n_bot = shape[0] // 2
    rows = [spec[:n_top]]
    if n_bot:
        rows.append(spec[full_shape[0] - n_bot:])
    return np.concatenate(rows, axis=0)[:, :shape[1] // 2 + 1]


def _double(arr, axis):
    """Exact 2x periodic linear upsampling of ``arr`` along ``axis``."""
    def sl(sub):
        return (sub, slice(None)) if axis == 0 else (slice(None), sub)

    shape = list(arr.shape)
    shape[axis] *= 2
    out = np.empty(shape, dtype=np.float32)
    lo, hi = arr[sl(slice(None, -1))], arr[sl(slice(1, None))]
    diff = hi - lo
    diff *= 0.25
    np.add(lo, diff, out=out[sl(slice(1, -1, 2))])
    diff *= 3
    np.add(lo, diff, out=out[sl(slice(2, None, 2))])
    first, last = arr[sl(slice(0, 1))], arr[sl(slice(-1, None))]
    out[sl(slice(0, 1))] = 0.75 * first + 0.25 * last
    out[sl(slice(-1, None))] = 0.75 * last + 0.25 * first
    return out


def _upsample(arr, shape):
    """
    Bilinear, periodic upsampling of ``arr`` onto ``shape``.

    ``arr.shape`` must lie on the pyramid of ``shape``, i.e. be
    ``_downsampled_shape(shape, 2**level)``. Upsampling is done by repeated
    exact doubling (cropping one pixel for odd sizes), which is much faster
    than ``ndimage.zoom`` or sparse interpolation for large arrays.
    """
    for axis in range(2):
        chain = [shape[axis]]
        while chain[-1] > arr.shape[axis]:
            chain.append(-(-chain[-1] // 2))
        if chain[-1] != arr.shape[axis]:
            raise ValueError(f"Cannot upsample {arr.shape} onto {shape}")
        for n_out in reversed(chain[:-1]):
            arr = _double(arr, axis)
            if arr.shape[axis] != n_out:
                arr = np.ascontiguousarray(
                    arr[:n_out] if axis == 0 else arr[:, :n_out])
    return arr


def _pyramid_factor(shape, min_size):
    """Largest power of 2 that keeps the short side >= ``min_size``."""
    return 2 ** max(0, int(np.log2(max(1, min(shape) // min_size))))


def _downsampled_shape(shape, factor):
    return tuple(int(np.ceil(n / factor)) for n in shape)


def _power_law_field(white, kk, shape, slope, k_outer, inner_scale):
    """Gaussian random field with P(k) ~ (k^2 + k_outer^2)^(-slope/2)."""
    amp = (kk**2 + k_outer**2) ** (-slope / 4)
    if inner_scale > 0:
        amp *= np.exp(-2 * (np.pi * inner_scale * kk) ** 2)
    amp[0, 0] = 0
    return _normalise(_irfft2(white * amp, shape))


def _smooth_field(shape, rng, slope=6.0, min_size=64):
    """
    Very smooth unit-variance random field (steep spectrum).

    It has no small-scale power, so it is computed on a coarse grid
    (``min_size`` pixels on the short side) and upsampled.
    """
    factor = _pyramid_factor(shape, min_size)
    shape_lo = _downsampled_shape(shape, factor)
    kk = _kgrid(shape_lo, shape)
    field = _power_law_field(_white_spectrum(shape_lo, rng), kk, shape_lo,
                             slope, 1 / min(shape), 0)
    return _upsample(field, shape)


def _ridged_field(white, kk, shape, slope, k_outer, inner_scale, sharpness):
    """
    Multi-octave ridged noise built from octave bands of ``white``.

    Each band is normalised and turned into ridges with exp(-p|g|). Octaves
    are weighted like a power-law field of the same slope. Each octave is
    also modulated by the ridges of the coarser octaves (heterogeneous
    multifractal), so fine filaments cluster along large-scale ones.

    Coarse octaves are synthesised on downsampled grids (at least ~8 pixels
    per wavelength) and upsampled as the loop moves to finer octaves.
    """
    k_min = max(k_outer, 1 / min(shape))
    k_max = 0.5
    if inner_scale > 0:
        k_max = min(k_max, 1 / (np.pi * inner_scale))
    n_oct = max(1, int(np.floor(np.log2(k_max / k_min))))
    max_level = np.log2(max(1, min(shape) // 16))

    with np.errstate(divide="ignore"):
        logk_full = np.log2(kk)
    total = mask = None
    for i in range(n_oct):
        k_c = k_min * 2 ** (i + 0.5)
        level = int(np.clip(np.floor(np.log2(0.125 / k_c)), 0, max_level))
        shape_lo = _downsampled_shape(shape, 2 ** level)
        if total is None:
            total = np.zeros(shape_lo, dtype=np.float32)
            mask = np.ones(shape_lo, dtype=np.float32)
        elif total.shape != shape_lo:
            total = _upsample(total, shape_lo)
            mask = _upsample(mask, shape_lo)

        # Smooth log-Gaussian band window, ~1 octave wide, avoids ringing
        logk = _crop_spectrum(logk_full, shape, shape_lo)
        window = np.exp(-2 * (logk - np.float32(np.log2(k_c))) ** 2)
        spec = _crop_spectrum(white, shape, shape_lo)
        ridge = _normalise(_irfft2(spec * window, shape_lo))
        np.abs(ridge, out=ridge)
        ridge *= -sharpness
        np.exp(ridge, out=ridge)

        ridge *= mask
        total += np.float32(k_c ** ((2 - slope) / 2)) * ridge
        if i < n_oct - 1:
            ridge *= 2.5
            mask = np.minimum(ridge, 1, out=ridge)

    return _normalise(_upsample(total, shape))


def _warp(field, rng, shape, strength):
    """Displace ``field`` by a smooth random vector field."""
    amp = min(shape) * strength
    coords = np.indices(shape, dtype=np.float32)
    for axis in range(2):
        coords[axis] += amp * _smooth_field(shape, rng)

    # map_coordinates is single-threaded but releases the GIL, so
    # interpolate horizontal strips in parallel
    out = np.empty(shape, dtype=np.float32)
    n_threads = min(8, os.cpu_count() or 1)
    edges = np.linspace(0, shape[0], 2 * n_threads + 1).astype(int)
    strips = [slice(lo, hi) for lo, hi in zip(edges[:-1], edges[1:]) if hi > lo]

    def interp(strip):
        ndimage.map_coordinates(field, coords[:, strip], output=out[strip],
                                order=1, mode="grid-wrap")

    with ThreadPoolExecutor(n_threads) as pool:
        list(pool.map(interp, strips))
    return out


def _envelope(kind, full_shape, rng, radius=0.35, width=0.05, ellipticity=0.2,
              angle=0.0, center=(0.5, 0.5), roughness=0.15, interior=0.1,
              asymmetry=0.5):
    """
    Large-scale brightness envelope.

    All lengths are fractions of ``min(shape)``, and ``center`` is
    (y, x) as fractions of the image. ``angle`` is in degrees.
    """
    # The envelope is smooth, so evaluate it on a coarse grid and upsample
    factor = _pyramid_factor(full_shape, 512)
    shape = _downsampled_shape(full_shape, factor)
    ny, nx = shape
    size = min(shape)
    yy, xx = np.indices(shape, dtype=np.float32)
    yy = (yy - center[0] * ny) / size
    xx = (xx - center[1] * nx) / size
    ang = np.deg2rad(angle)
    xr = xx * np.cos(ang) + yy * np.sin(ang)
    yr = (-xx * np.sin(ang) + yy * np.cos(ang)) / (1 - ellipticity)
    rr = np.hypot(xr, yr)

    if kind == "blob":
        env = np.exp(-2 * (rr / radius) ** 2)  # radius = 2 sigma
        return _upsample(env, full_shape)

    if kind == "shell":
        bumps = _smooth_field(shape, rng)
        r_shell = radius * (1 + roughness * bumps)
        env = np.exp(-0.5 * ((rr - r_shell) / width) ** 2)
        env += interior * (rr < r_shell)
        # Brighter on one side, like a shell expanding into a density gradient
        cos_theta = np.divide(xr, rr, out=np.zeros_like(rr), where=rr > 0)
        env *= 1 + asymmetry * cos_theta
        return _upsample(env, full_shape)

    raise ValueError(f"Unknown envelope '{kind}', use 'shell', 'blob', "
                     "an array, or None.")


def nebulosity(shape=(1024, 1024), *, power_slope=2.7, outer_scale=0.5,
               inner_scale=0.0, filament_strength=0.5, filament_sharpness=3.0,
               warp_strength=0.03, contrast=0.8, envelope=None,
               envelope_kwargs=None, seed=None, dtype=np.float32):
    """
    Generate a semi-realistic 2D nebulosity image.

    Lognormal turbulence, ridged filaments and domain warping are combined
    using only FFTs and one interpolation. A 2048 x 2048 image takes about
    one second on a laptop. Without an envelope the result is periodic, so
    it tiles seamlessly.

    Parameters
    ----------
    shape : tuple of int
        (ny, nx) of the output array.
    power_slope : float
        Slope beta of the power spectrum, P(k) ~ k^-beta. Measured ISM
        values are about 2.6-3.0. Larger values give smoother images.
    outer_scale : float
        Largest structure size as a fraction of ``min(shape)``. The
        spectrum flattens below the corresponding wavenumber.
    inner_scale : float
        Gaussian smoothing scale in pixels. Use 0 for no smoothing.
    filament_strength : float
        Mixing fraction [0, 1] of the ridged filament component.
    filament_sharpness : float
        Ridge exponent. Larger values give thinner filaments.
    warp_strength : float
        Domain-warp displacement amplitude as a fraction of ``min(shape)``.
        Use 0 to disable.
    contrast : float
        Standard deviation of the log-intensity, i.e. the lognormal sigma.
        Measured values for NGC 6888 are about 0.4-0.8. Above about 1.5 a
        few bright knots dominate the image.
    envelope : {None, "shell", "blob"} or array_like
        Large-scale brightness profile multiplied into the image. An array
        must have the same ``shape``.
    envelope_kwargs : dict, optional
        Options for the built-in envelopes: ``radius``, ``width``,
        ``ellipticity``, ``angle`` (deg), ``center`` ((y, x) fractions),
        ``roughness``, ``interior``, ``asymmetry``. Lengths are fractions
        of ``min(shape)``. For ``"blob"``, ``radius`` is the 2-sigma
        radius of the Gaussian.
    seed : int, np.random.Generator, optional
        Random seed. Defaults to ``__config__["!random.seed"]``.
    dtype : numpy dtype
        Output dtype.

    Returns
    -------
    image : np.ndarray
        Non-negative array with maximum 1.

    Examples
    --------
    ::

        from astropy import units as u
        from synphot import SourceSpectrum, ConstFlux1D
        from scopesim_templates.misc import nebulosity, source_from_array

        img = nebulosity((512, 512), envelope="shell", seed=42)
        sed = SourceSpectrum(ConstFlux1D, amplitude=1)
        src = source_from_array(img, sed=sed, pixel_scale=0.1*u.arcsec,
                                amplitude=15*u.ABmag, filter_curve="Ks")
    """
    shape = tuple(int(n) for n in shape)
    if len(shape) != 2:
        raise ValueError(f"shape must be 2D, got {shape}")
    if not 0 <= filament_strength <= 1:
        raise ValueError("filament_strength must be in [0, 1]")
    if seed is None:
        seed = __config__["!random.seed"]
    rng = np.random.default_rng(seed)

    kk = _kgrid(shape)
    k_outer = 1 / (outer_scale * min(shape))
    white = _white_spectrum(shape, rng)

    log_img = _power_law_field(white, kk, shape, power_slope, k_outer,
                               inner_scale)
    if filament_strength > 0:
        ridged = _ridged_field(white, kk, shape, power_slope, k_outer,
                               inner_scale, filament_sharpness)
        log_img = ((1 - filament_strength) * log_img
                   + filament_strength * ridged)
    if warp_strength > 0:
        log_img = _warp(log_img, rng, shape, warp_strength)

    image = np.exp(contrast * _normalise(log_img))

    if envelope is not None:
        if isinstance(envelope, str):
            env = _envelope(envelope, shape, rng,
                            **(envelope_kwargs or {}))
        else:
            env = np.asarray(envelope)
            if env.shape != shape:
                raise ValueError(f"envelope shape {env.shape} does not match "
                                 f"{shape}")
        image *= env

    image = np.clip(image, 0, None)
    image /= image.max()
    return image.astype(dtype, copy=False)
