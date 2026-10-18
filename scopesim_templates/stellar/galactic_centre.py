# -*- coding: utf-8 -*-
"""Galactic Centre star-field template (Gillessen+2017 orbits)."""

import warnings

import numpy as np
from astropy import units as u
from astropy.table import Table

import pyckles
from spextra import Spextrum

from ..rc import Source
from ..utils.general_utils import add_function_call_str
from . import _gillessen
from . import stars_utils as su
from ._orbits import RA_SGR_A_STAR, DEC_SGR_A_STAR
from .gc_background import galactic_centre_background


__all__ = ["galactic_centre"]


@add_function_call_str
def galactic_centre(time,
                    filter_name="Paranal/NACO.Ks",
                    early_type_spec="B0V",
                    late_type_spec="M8III",
                    background_mass=0.0,
                    seed=None,
                    background_kwargs=None):
    """
    Source of the Sgr A* star field at a user-specified epoch.

    Stars and orbital elements from Gillessen et al. 2017
    (CDS J/ApJ/837/30, vendored). Per-star spectra are RV-shifted to the
    line-of-sight velocity at ``time`` and scaled to each star's catalogue
    K-band magnitude.

    Parameters
    ----------
    time : astropy.time.Time, datetime, ISO string, or float MJD
        Epoch at which to evaluate the orbits.
    filter_name : str
        SVO filter ID used to scale the spectra to ``Kmag``. Defaults to
        ``"Paranal/NACO.Ks"`` to match the photometric system of the
        Gillessen catalogue.
    early_type_spec, late_type_spec : str
        Spectral types assigned to early-type (``SpT == "e"``) and
        late-type (``SpT == "l"``) stars respectively.
    background_mass : float or str
        [Msun] If > 0, add a hypothetical unseen main-sequence population
        of this total mass on random orbits about Sgr A*, built by
        :func:`~scopesim_templates.stellar.galactic_centre_background`.
        Default 0 (S-stars only). ``"gravity2022"`` selects the
        3000 Msun upper limit of GRAVITY Collaboration 2022 (A&A 657, L12).
    seed : int, np.random.Generator, or None
        Seed for the background population; ignored if
        ``background_mass == 0``. Re-use an integer seed at different
        ``time`` values to follow the same background stars through time.
    background_kwargs : dict, optional
        Further keyword arguments for
        :func:`~scopesim_templates.stellar.galactic_centre_background`,
        e.g. ``{"r_max": 1.0, "mag_limit": 24}``.

    Returns
    -------
    scopesim.Source
        With ``background_mass > 0`` the Source has two fields: the
        S-stars (``fields[0]``) and the background population
        (``fields[1]``).

    Examples
    --------
    >>> from astropy.time import Time
    >>> from scopesim_templates.stellar import galactic_centre
    >>> src = galactic_centre(Time("2024-01-01"))

    With the GRAVITY 2022 upper limit of unseen stars added:

    >>> src = galactic_centre(Time("2031-10-22"), background_mass=3000,
    ...                       seed=42, background_kwargs={"mag_limit": 24})
    """
    tbl = _gillessen.stars_at_time(time)
    n_stars = len(tbl)

    # Map SpT codes to user-supplied spectral types. The Gillessen catalogue
    # has empty SpT cells for a small number of stars (e.g. S39, S55) that
    # are otherwise valid; warn and treat those as early-type.
    spt_map = {"e": early_type_spec, "l": late_type_spec}
    spec_types = []
    for code in tbl["SpT"]:
        key = str(code).strip().lower()
        spt = spt_map.get(key)
        if spt is None:
            warnings.warn(
                f"galactic_centre: missing or unrecognised SpT {code!r}; "
                f"defaulting to early-type {early_type_spec!r}",
                UserWarning,
                stacklevel=2,
            )
            spt = early_type_spec
        spec_types.append(spt)

    # Base spectrum per unique user-supplied spectral type, scaled to 0 mag.
    unique_user_types = [early_type_spec, late_type_spec]
    pickles_lib = pyckles.SpectralLibrary("pickles", return_style="synphot")
    cat_types = su.nearest_spec_type(unique_user_types, pickles_lib.table)
    base = {
        orig: Spextrum(modelclass=pickles_lib[cat]).scale_to_magnitude(
            filter_curve=filter_name, amplitude=0 * u.mag)
        for orig, cat in zip(unique_user_types, cat_types)
    }

    # Per-star spectra: each gets its own RV-shifted copy of the base.
    # Pass rv as a Quantity: Spextrum.redshift(vel=...) assumes m/s for
    # bare floats, which would silently under-apply the km/s shifts by 1000.
    rv = tbl["RV"].to(u.km / u.s)
    spectra = [base[spec_types[i]].redshift(vel=rv[i])
               for i in range(n_stars)]

    # Weights bake in the K-band magnitude (base is at 0 mag).
    weight = 10 ** (-0.4 * tbl["Kmag"].to(u.mag).value)

    # Identity ref: each star points at its own (unique) spectrum.
    ref = np.arange(n_stars, dtype=int)

    field = Table(
        names=["x", "y", "ref", "weight", "rv", "spec_types"],
        data=[tbl["d_ra"].to(u.arcsec),
              tbl["d_dec"].to(u.arcsec),
              ref,
              weight,
              rv,
              spec_types],
    )

    src = Source(spectra=spectra, table=field)
    src.meta.update({
        "description": "Galactic Centre (Gillessen+2017)",
        "epoch": str(time),
        "ra": RA_SGR_A_STAR,
        "dec": DEC_SGR_A_STAR,
        "filter_name": filter_name,
    })

    if isinstance(background_mass, str) or background_mass > 0:
        background = galactic_centre_background(
            time, total_mass=background_mass, filter_name=filter_name,
            seed=seed, **(background_kwargs or {}))
        meta = dict(src.meta)
        meta.update({
            "background_mass": background.meta["total_mass"],
            "background_n_stars": len(background.fields[0].field),
            "seed": background.meta["seed"],
        })
        src = src + background
        src.meta.clear()
        src.meta.update(meta)
    return src
