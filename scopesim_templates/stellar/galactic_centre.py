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


__all__ = ["galactic_centre"]


# Sgr A* (J2000)
RA_SGR_A_STAR = 266.41684  # deg
DEC_SGR_A_STAR = -29.00781  # deg


@add_function_call_str
def galactic_centre(time,
                    filter_name="Paranal/NACO.Ks",
                    early_type_spec="B0V",
                    late_type_spec="M8III"):
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

    Returns
    -------
    scopesim.Source

    Examples
    --------
    >>> from astropy.time import Time
    >>> from scopesim_templates.stellar import galactic_centre
    >>> src = galactic_centre(Time("2024-01-01"))
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
        "object": "Galactic Centre (Gillessen+2017)",
        "epoch": str(time),
        "ra": RA_SGR_A_STAR,
        "dec": DEC_SGR_A_STAR,
        "filter_name": filter_name,
    })
    return src
