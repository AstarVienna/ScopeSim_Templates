import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.utils.decorators import deprecated_renamed_argument

import pyckles

from .. import utils
from .. import rc
from ..rc import ter_curve_utils as tcu


def star_field(n, mmin, mmax, width, height=None, filter_name="V", **kwargs):
    """
    Creates a super basic field of stars with random positions and brightnesses

    Parameters
    ----------
    n : int
        number of stars
    mmin, mmax : astropy.Quantity, float
        [mag, Jy] minimum and maximum flux amplitudes of the population
        in u.mag (u.ABmag) or u.Janksy where u is for astropy.units.
        If mmin and mmax are floats, Vega magnitudes are assumed
    width, height : float
        [arcsec] width of region to put stars in. if height=None, height=width
    filter_name : str
        For scaling the stars. Use either common names or Spanish-VO identifiers

    kwargs
    ------
    x, y : lists, arrays
        [arcsec] The positions of the stars can be overridden by specifying the
        coordinates. The lists must contain N values

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

    if rc.__config__["!random.seed"] is not None:
        np.random.seed(rc.__config__["!random.seed"])

    if not isinstance(mmin, u.Quantity) and not isinstance(mmax, u.Quantity):
        mmin, mmax = u.Quantity(mmin, u.mag), u.Quantity(mmax, u.mag)

    if "x" in kwargs and "y" in kwargs:
        x, y = kwargs["x"], kwargs["y"]
    else:
        rands = np.random.random(size=(2, n)) - 0.5
        x, y = width * rands[0], height * rands[1]

    # amplitudes = np.random.random(size=n) * (mmax - mmin) + mmin
    amplitudes = np.linspace(mmin, mmax, n)
    spec_types = ["A0V"] * n

    src = stars(filter_name=filter_name, amplitudes=amplitudes,
                spec_types=spec_types, x=x, y=y)
    src.meta["scaling_unit"] = mmin.unit

    return src


def star_grid(n, mmin, mmax, filter_name="V", separation=1):
    """
    Creates a square grid of A0V stars at equal magnitude intervals

    Parameters
    ----------
    n : int
        the number of stars in the grid
    mag_min, mag_max : float
        [vega mag] the minimum (brightest) and maximum (faintest) magnitudes for
        stars in the grid
    filter_name : str
        any filter that is in the SimCADO package directory.
        See ``scopesim.optics.get_filter_set()``
    separation : float, optional
        [arcsec] an average speration between the stars in the grid can be
        specified. Default is 1 arcsec

    Returns
    -------
    src : ``scopesim.Source``

    """
    side_len = int(np.ceil(np.sqrt(n)))
    x = separation * (np.arange(n) % side_len - (side_len - 1) / 2)
    y = separation * (np.arange(n) // side_len - (side_len - 1) / 2)
    src = star_field(n, mmin, mmax, side_len, filter_name=filter_name, x=x, y=y)

    return src


@deprecated_renamed_argument('mags', 'amplitudes', '0.1')
def stars(filter_name, amplitudes, spec_types, x, y):
    """
    Creates a scopesim.Source object for a list of stars with given amplitudes

    .. note:: If amplitudes have no units, vega magnitudes are assumed

    Parameters
    ----------
    filter_name : str
        For scaling the stars. Use either common names or Spanish-VO identifiers
    amplitudes : list of Quanitity, float
        [mag, Jy] amplitudes for the list of stars. Acceptable astropy.units:
        [u.mag, u.ABmag, u.Janksy]. If no units are given, Vega magnitudes are
        assumed
    spec_types : list of strings
        the spectral type(s) of the stars, e.g. "A0V", "G5III"
    x, y : arrays of float
        [arcsec] x and y coordinates of the stars on the focal plane

    Returns
    -------
    src : ``scopesim.Source``


    Examples
    --------

    Create a ``Source`` object for a random group of stars::

        >>> import numpy as np
        >>> import astropy.units as u
        >>> from scopesim_templates.basic.stars import stars
        >>>
        >>> n = 100
        >>> spec_types = ["A0V", "G2V", "K0III", "M5III", "O8I"]
        >>> ids = np.random.randint(0, 5, size=n)
        >>> star_list = [spec_types[i] for i in ids]
        >>> mags = np.random.normal(loc=20, scale=3, size=n) * u.mag
        >>> x, y = np.random.random(size=(2, n))
        >>>
        >>> src = stars("Ks", mags, spec_types, x, y)

    **All positions are in arcsec.**

    The final positions table is kept in the ``<Source>.fields`` attribute::

        >>> src.fields[0]

    Each star in this table has an associated spectrum kept in the
    ``<Source>.spectra`` attribute. These stars are connected to the spectra in
    the list by the "ref" column in the ``.fields`` table::

        >>> src.spectra

    The stars can be scaled in units of u.mag, u.ABmag or u.Jansky. Any filter
    listed on the spanish VO filter profile service can be used for the scaling
    (``http://svo2.cab.inta-csic.es/theory/fps/``). SVO filter names need to be
    in the following format ``observatory/instrument.filter``::

        >>> amplitudes = np.linspace(1, 3631, n) * u.Jansky
        >>> filter_name = "Paranal/HAWKI.Ks"
        >>> stars(filter_name, amplitudes, spec_types, x=x, y=y)

    """
    if not isinstance(spec_types, (list, tuple, np.ndarray)):
        spec_types = [spec_types]

    if not isinstance(amplitudes, u.Quantity):
        amplitudes = u.Quantity(amplitudes, u.mag, copy=False)

    if not isinstance(x, u.Quantity):
        x = u.Quantity(x, u.arcsec, copy=False)
    if not isinstance(y, u.Quantity):
        y = u.Quantity(y, u.arcsec, copy=False)

    pickles_lib = pyckles.SpectralLibrary("pickles", return_style="synphot")
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

    spectra = [tcu.scale_spectrum(pickles_lib[spt], filter_name, zero)
               for spt in zip(cat_spec_types)]

    # get the references to the unique stellar types
    ref_dict = {spt: ii for ii, spt in enumerate(unique_types)}
    ref = np.array([ref_dict[i] for i in spec_types])

    src = rc.Source(spectra=spectra, x=x, y=y, ref=ref, weight=weight)

    return src


def star(filter_name, amplitude, spec_type="A0V", x=0, y=0):
    src = stars(filter_name, [amplitude], [spec_type], [x], [y])
    return src


star.__doc__ = stars.__doc__
