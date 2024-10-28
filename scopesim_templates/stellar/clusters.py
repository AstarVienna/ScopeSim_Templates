"""Currently only contains a simple zero-age open cluster."""

import numpy as np
import pyckles
from astropy import units as u
from astropy.table import Table
from astropy.utils import deprecated_renamed_argument

from spextra import Spextrum

from ..rc import __config__, Source
from ..utils.general_utils import function_call_str, RA0, DEC0
from . import imf, cluster_utils as cu


__all__ = ["cluster"]


@deprecated_renamed_argument("half_light_radius", "core_radius", "0.1")
def cluster(mass=1E3, distance=50000, core_radius=1, ra=RA0, dec=DEC0,
            **kwargs):
    """
    Generate a source object for a young cluster.

    The cluster distribution follows a gaussian profile with the
    ``core_radius`` corresponding to the HWHM of the distribution. The
    choice of stars follows a Kroupa IMF, with no evolved stars in the mix.
    Ergo this is more suitable for a young cluster than an evolved custer


    Parameters
    ----------
    mass : float
        [Msun] Mass of the cluster (not number of stars). Max = 1E5 Msun
    distance : float
        [pc] distance to the cluster
    core_radius : float
        [pc] half light radius of the cluster
    ra : float, str
        RA of the source
    dec : float, str
        DEC of the source

    Additional parameters
    ---------------------
    tidal_radius : float
        [pc] Not yet implemented, for later once there is a King profile
    multiplicity : Unknown
        A multiplicity object from the imf.py package. Not sure what it does.
    seed: float
        For setting the random seed for both masses and positions

    Returns
    -------
    src: scopesim.Source


    Examples
    --------
    Create a ``Source`` object for a young open cluster with half light radius
    of around 0.2 pc at the extragalactic centre and 1000 solar masses worth
    of stars:

        >>> from scopesim_templates.stellar.clusters import cluster
        >>> src = cluster(mass=1000, distance=8500, core_radius=0.2, seed=9001)

    """
    params = {"mass": mass,
              "distance": distance,
              "core_radius": core_radius,
              "tidal_radius": None,
              "multiplicity_object": None,
              "seed": __config__["!random.seed"],
              "ra": ra,
              "dec": dec}
    params.update(kwargs)
    params["function_call"] = function_call_str(cluster, kwargs=params)
    params["object"] = "cluster"

    # 1. sample masses from an IMF
    kroupa = imf.Kroupa_2001(params["multiplicity_object"])
    masses, _, _, _ = kroupa.generate_cluster(mass, seed=params["seed"])
    masses.clip(max=250, out=masses)

    # 2. get spec_types for masses
    spec_types = cu.mass2spt(masses)
    spec_types = cu.closest_pickles(spec_types)

    # 3. get spectra from pyckles
    pickles = pyckles.SpectralLibrary("pickles", return_style="synphot")
    unique_spts = set(spec_types)
    spectra = [pickles[spt] for spt in unique_spts]

    # 4. scale all spectra to V=0
    # spectra = [tcu.scale_spectrum(spec, "V", 0*u.mag) for spec in spectra]
    spectra = [Spextrum(modelclass=spec).scale_to_magnitude(filter_curve="V",
                                                            amplitude=0*u.mag)
               for spec in spectra]

    # 5. make ref list for spec_types from spectra
    ref = [np.where([spt == u_spt for u_spt in unique_spts])[0][0]
           for spt in spec_types]

    # 6. make weight list from Mv + dist_mod(distance)
    Mvs = np.array(cu.mass2absmag(masses))
    dist_mod = 5 * np.log10(distance) - 5
    weight = 10 ** (-0.4 * (Mvs + dist_mod))

    # 7. make x,y from half_light_radius and distance
    rad2arcsec = 206264.80624709636
    fwhm = 2 * core_radius / distance * rad2arcsec
    x, y = cu.gaussian_distribution(len(masses), fwhm=fwhm,
                                    seed=params["seed"])
    x = x << u.arcsec
    y = y << u.arcsec

    # 8. make table with (x,y,ref,weight)
    tbl = Table(names=["x", "y", "ref", "weight", "masses", "spec_types"],
                data=[x, y, ref, weight, masses, spec_types],
                units=[u.arcsec, u.arcsec, None, None, u.solMass, None])

    # 9. make Source with table, spectra
    src = Source(table=tbl, spectra=spectra)
    src.meta.update(params)

    return src
