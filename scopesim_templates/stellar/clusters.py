import numpy as np
import pyckles
from astropy import units as u
from astropy.table import Table
from astropy.utils import deprecated_renamed_argument
from scopesim_templates import rc
from scopesim_templates.rc import ter_curve_utils as tcu
from scopesim_templates.stellar import imf, cluster_utils as cu
from scopesim_templates.utils.general_utils import function_call_str


@deprecated_renamed_argument('half_light_radius', 'core_radius', '0.1')
def cluster(mass=1E3, distance=50000, core_radius=1, **kwargs):
    """
    Generate a source object for a young cluster

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
    of around 0.2 pc at the extragalactic centre and 1000 solar masses worth of stars:

        >>> from scopesim_templates.basic.stars import cluster
        >>> src = cluster(mass=1000, distance=8500, core_radius=0.2, seed=9001)

    """
    params = {"mass": mass,
              "distance": distance,
              "core_radius": core_radius,
              "tidal_radius": None,
              "multiplicity_object": None,
              "seed": rc.__config__["!random.seed"]}
    params.update(kwargs)
    params["function_call"] = function_call_str(cluster, params)
    params["object"] = "cluster"

    # 1. sample masses from an IMF
    kroupa = imf.Kroupa_2001(params["multiplicity_object"])
    results = kroupa.generate_cluster(mass, seed=params["seed"])
    masses, multi_flag, secondary_masses, system_masses = results
    masses[masses > 250] = 250

    # 2. get spec_types for masses
    spec_types = cu.mass2spt(masses)
    spec_types = cu.closest_pickles(spec_types)

    # 3. get spectra from pyckles
    pickles = pyckles.SpectralLibrary("pickles", return_style="synphot")
    unique_spts = np.unique(spec_types)
    spectra = [pickles[spt] for spt in unique_spts]

    # 4. scale all spectra to V=0
    spectra = [tcu.scale_spectrum(spec, "V", 0*u.mag) for spec in spectra]

    # 5. make ref list for spec_types from spectra
    ref = [np.where([spt == u_spt for u_spt in unique_spts])[0][0]
           for spt in spec_types]

    # 6. make weight list from Mv + dist_mod(distance)
    Mvs = np.array(cu.mass2Mv(masses))
    dist_mod = 5 * np.log10(distance) - 5
    weight = 10 ** (-0.4 * (Mvs + dist_mod))

    # 7. make x,y from half_light_radius and distance
    rad2arcsec = 206264.80624709636
    fwhm = 2 * core_radius / distance * rad2arcsec
    x, y = cu.gaussian_distribution(len(masses), fwhm=fwhm, seed=params["seed"])

    # 8. make table with (x,y,ref,weight)
    tbl = Table(names=["x", "y", "ref", "weight", "masses", "spec_types"],
                data= [ x,   y,   ref,   weight,   masses,   spec_types ])

    # 9. make Source with table, spectra
    src = rc.Source(table=tbl, spectra=spectra)
    src.meta.update(params)

    return src