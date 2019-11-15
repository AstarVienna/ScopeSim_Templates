import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.utils.decorators import deprecated_renamed_argument

import pyckles

from .. import utils
from .. import rc


def star_field(n, mag_min, mag_max, width, height=None,
               photometric_system="vega"):
    """
    Creates a super basic field of stars with random positions and brightnesses

    Parameters
    ----------
    n : int
        number of stars

    mag_min, mag_max : float
        [mag] minimum and maximum magnitudes of the population

    width : float
        [arcsec] width of region to put stars in

    height : float, optional
        [arcsec] if None, then height=width

    photometric_system : str, optional
        [vega, AB]


    Returns
    -------
    stars : scopesim.Source object
        A Source object with a field of stars that can be fed into the method:
        ``<OpticalTrain>.observe()``

    See Also
    --------
    ``<OpticalTrain>.observe``
    ``<OpticalTrain>.readout``

    """
    if height is None:
        height = width

    if photometric_system.lower() == "ab":
        spec = utils.ab_spectrum()
    else:
        spec = utils.vega_spectrum()

    if rc.__config__["!SIM.random.seed"] is not None:
        np.random.seed(rc.__config__["!SIM.random.seed"])

    rands = np.random.random(size=(2, n)) - 0.5
    x = width * rands[0]
    y = height * rands[1]
    mags = np.random.random(size=n) * (mag_max - mag_min) + mag_min
    w = 10**(-0.4 * mags)
    ref = np.zeros(n, dtype=int)

    tbl = Table(data=[x, y, w, ref, mags],
                names=["x", "y", "weight", "ref", "mag"])
    tbl.meta["photometric_system"] = photometric_system
    stars = rc.Source(spectra=spec, table=tbl)

    return stars


@deprecated_renamed_argument('mags', 'amplitudes', '0.1')
def stars(spec_types, amplitudes, filter_name, x, y):
    """
    If amplitude is a float, it is assumed to by a vega magnitude

    Parameters
    ----------
    spec_types : list of str
    amplitude : list of Quanitity, float
    filter_name : str
    x, y : list of float

    Returns
    -------

    """
    if not isinstance(spec_types, (list, tuple, np.ndarray)):
        spec_types = [spec_types]

    if not isinstance(amplitudes, u.Quantity):
        amplitudes = u.Quantity(amplitudes, u.mag, copy=False)

    if not isinstance(x, u.Quantity):
        x = u.Quantity(x, u.arcsec, copy=False)
    if not isinstance(y, u.Quantity):
        y = u.Quantity(y, u.arcsec, copy=False)

    pickles_lib = pyckles.load_catalog("pickles")
    unique_types = np.unique(spec_types)
    cat_spec_types = utils.nearest_spec_type(unique_types,
                                             pickles_lib.table)

    # scale the spectra and get the weights
    if amplitudes.unit in [u.mag, u.ABmag, u.STmag]:
        zero = 0 * amplitudes.unit
        weight = 10**(-0.4*amplitudes.value)
    else:
        zero = 1 * amplitudes.unit
        weight = amplitudes.value


    # ..todo FIND TER_CUREVE_UTILS


    spectra = [rc.scale_spectrum(pickles_lib[spt], zero, filter_name)
               for spt in zip(cat_spec_types)]

    # get the references to the unique stellar types
    ref_dict = {spt: ii for ii, spt in enumerate(unique_types)}
    ref = np.array([ref_dict[i] for i in spec_types])

    src = rc.Source(spectra=spectra, x=x, y=y, ref=ref, weight=weight)

    return src





