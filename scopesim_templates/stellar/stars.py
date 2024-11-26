import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.utils.decorators import deprecated_renamed_argument

import pyckles
from spextra import Spextrum

from ..rc import __config__, Source
from ..utils.general_utils import RA0, DEC0, add_function_call_str,\
    function_call_str
from . import stars_utils as su


__all__ = ["star",
           "stars",
           "star_field",
           "star_grid",
           ]


def star_field(n, mmin, mmax, width, height=None, filter_name="V",
               ra=RA0, dec=DEC0, **kwargs):
    """
    Create a super basic field of stars with random positions and brightnesses.

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
        For scaling the stars. Use either common names or Spanish-VO
        identifiers.
    ra : float or str
        RA of the center of the field (not used at the moment).
    dec : float or str
        DEC of the center of the field (not used at the moment).

    Additional parameters
    ---------------------
    x, y : lists, arrays
        [arcsec] The positions of the stars can be overridden by specifying the
        coordinates. The lists must contain N values

    Returns
    -------
    stars : scopesim.Source object
        A Source object with a field of stars that can be fed into the method:
        ``<OpticalTrain>.observe()``

    """
    params = {"n": n,
              "mmin": mmin,
              "mmax": mmax,
              "width": width,
              "height": height,
              "filter_name": "V",
              "seed": __config__["!random.seed"],
              "ra": ra,
              "dec": dec}
    params.update(kwargs)
    params["function_call"] = function_call_str(star_field, kwargs=params)
    params["object"] = "star field"

    if height is None:
        height = width

    if isinstance(params["seed"], int):
        np.random.seed(params["seed"])

    if not isinstance(mmin, u.Quantity) and not isinstance(mmax, u.Quantity):
        mmin, mmax = u.Quantity(mmin, u.mag), u.Quantity(mmax, u.mag)

    if "x" in kwargs and "y" in kwargs:
        x, y = kwargs.pop("x"), kwargs.pop("y")
    else:
        rands = np.random.random(size=(2, n)) - 0.5
        x, y = width * rands[0], height * rands[1]

    # amplitudes = np.random.random(size=n) * (mmax - mmin) + mmin
    amplitudes = np.linspace(mmin, mmax, n)
    spec_types = ["A0V"] * n

    src = stars(filter_name=filter_name, amplitudes=amplitudes,
                spec_types=spec_types, x=x, y=y, ra=ra, dec=dec, **kwargs)
    src.meta["scaling_unit"] = mmin.unit
    src.meta.update(params)

    return src


@add_function_call_str
def star_grid(n, mmin, mmax, filter_name="V", separation=1, ra=RA0, dec=DEC0):
    """
    Create a square grid of A0V stars at equal magnitude intervals.

    Parameters
    ----------
    n : int
        the number of stars in the grid
    mmin : float
        [vega mag] the minimum (brightest) magntude of stars in the grid
    mmax : float
        [vega mag] maximum (faintest) magnitudes for stars in the grid
    filter_name : str
        any filter that is in the SimCADO package directory.
        See ``scopesim.optics.get_filter_set()``
    separation : float, optional
        [arcsec] an average speration between the stars in the grid can be
        specified. Default is 1 arcsec
    ra : float or str
        RA of the center of the field (not used at the moment)
    dec : float or str
        DEC of the center of the field (not used at the moment)

    Returns
    -------
    src : ``scopesim.Source``

    """
    side_len = int(np.ceil(np.sqrt(n)))
    x = separation * (np.arange(n) % side_len - (side_len - 1) / 2)
    y = separation * (np.arange(n) // side_len - (side_len - 1) / 2)

    src = star_field(n, mmin, mmax, side_len, filter_name=filter_name,
                     x=x, y=y, ra=ra, dec=dec)
    return src


@deprecated_renamed_argument("mags", "amplitudes", "0.1")
@add_function_call_str
def stars(filter_name, amplitudes, spec_types, x, y, library="pyckles",
          ra=RA0, dec=DEC0):
    """
    Create a scopesim.Source object for a list of stars with given amplitudes.

    .. note:: If amplitudes have no units, vega magnitudes are assumed

    Parameters
    ----------
    filter_name : str
        For scaling the stars. Use either common names or Spanish-VO
        identifiers.
    amplitudes : list of Quanitity, float
        [mag, Jy] amplitudes for the list of stars. Acceptable astropy.units:
        [u.mag, u.ABmag, u.Janksy]. If no units are given, Vega magnitudes are
        assumed
    spec_types : list of strings
        the spectral type(s) of the stars, e.g. "A0V", "G5III"
    x, y : arrays of float
        [arcsec] x and y coordinates of the stars on the focal plane
    ra, dec : float
        coordinates of the center of the field
    library: str
        Library where the spectroscopic types are taken. By default are taken
        from the pickles library using the `pyckles` package. Other libraries
        available are kurucz, bosz/lr, bosz/mr, bosz/hr, etc for MIR coverage
        and different spectral resolutions. Please see the `spextra` package
        for more information

    Returns
    -------
    src : ``scopesim.Source``

    Examples
    --------
    Create a ``Source`` object for a random group of stars::

        >>> import numpy as np
        >>> import astropy.units as u
        >>> from scopesim_templates.stellar import stars
        >>>
        >>> n = 100
        >>> spec_types = ["A0V", "G2V", "K0III", "M5III", "O8I"] * (n // 5)
        >>> ids = np.random.randint(0, 5, size=n)
        >>> star_list = [spec_types[i] for i in ids]
        >>> mags = np.random.normal(loc=20, scale=3, size=n) * u.mag
        >>> x, y = np.random.random(size=(2, n))
        >>>
        >>> src = stars("Ks", mags, spec_types, x, y)

    .. note: All positions are in arcsec.

    The final positions table is kept in the ``<Source>.fields`` attribute::

        >>> my_pos_table = src.fields[0]

    Each star in this table has an associated spectrum kept in the
    ``<Source>.spectra`` attribute. These stars are connected to the spectra in
    the list by the "ref" column in the ``.fields`` table::

        >>> list(src.spectra.values())
        [SpextrumNone, SpextrumNone, SpextrumNone, SpextrumNone, SpextrumNone]

    The stars can be scaled in units of u.mag, u.ABmag or u.Jansky. Any filter
    listed on the spanish VO filter profile service can be used for the scaling
    (``http://svo2.cab.inta-csic.es/theory/fps/``). SVO filter names need to be
    in the following format ``observatory/instrument.filter``::

        >>> amplitudes = np.linspace(1, 3631, n) * u.Jansky
        >>> filter_name = "Paranal/HAWKI.Ks"
        >>> my_source = stars(filter_name, amplitudes, spec_types, x=x, y=y)

    """
    if not isinstance(spec_types, (list, tuple, np.ndarray)):
        spec_types = [spec_types]

    if not isinstance(amplitudes, u.Quantity):
        amplitudes = u.Quantity(amplitudes, u.mag, copy=False)

    if not isinstance(x, u.Quantity):
        x = u.Quantity(x, u.arcsec, copy=False)
    if not isinstance(y, u.Quantity):
        y = u.Quantity(y, u.arcsec, copy=False)

    unique_types = list(set(spec_types))
    if library == "pyckles":
        pickles_lib = pyckles.SpectralLibrary("pickles",
                                              return_style="synphot")
        cat_spec_types = su.nearest_spec_type(unique_types, pickles_lib.table)
        # print(cat_spec_types)
    # scale the spectra and get the weights
    if amplitudes.unit in [u.mag, u.ABmag, u.STmag]:
        zero = 0 * amplitudes.unit
        weight = 10**(-0.4*amplitudes.value)
    else:
        zero = 1 * amplitudes.unit
        weight = amplitudes.value

    if library == "pyckles":
        # spectra = [tcu.scale_spectrum(pickles_lib[spt], filter_name, zero)
        #            for spt in zip(cat_spec_types)]
        spectra = [Spextrum(modelclass=pickles_lib[spt]).scale_to_magnitude(
            filter_curve=filter_name, amplitude=zero)
                   for spt in cat_spec_types]

    else:
        spectra = [Spextrum(library + "/" + spec.lower()).scale_to_magnitude(
            amp, filter_curve=filter_name)
                   for spec, amp in zip(spec_types, amplitudes)]
        # TODO: actually use weights here instead of spectra scaling to avoid
        #       O(n) scaling for number of spectra...
        weight = np.ones(shape=amplitudes.shape)

    # get the references to the unique stellar types
    ref = np.array([unique_types.index(i) for i in spec_types])

    if len(ref) == 1:
        ref = ref[0] * np.ones(len(x), dtype=int)
        spec_types *= len(x)

    tbl = Table(names=["x", "y", "ref", "weight", "spec_types"],
                data=[x, y, ref, weight, spec_types])

    src = Source(spectra=spectra, table=tbl)
    return src


def star(filter_name, amplitude, spec_type="A0V", x=0, y=0, library="pyckles"):
    if isinstance(amplitude, u.Quantity) is False:
        amplitude = amplitude * u.mag

    src = stars(filter_name, [amplitude.value] * amplitude.unit,
                [spec_type], [x], [y], library=library)

    return src


star.__doc__ = stars.__doc__
