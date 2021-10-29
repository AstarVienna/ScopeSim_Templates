import numpy as np


def nearest_spec_type(spec, cat_tbl):
    """
    Finds the closest available spectral type in a catalogue

    Parameters
    ----------
    spec : str
        Stellar spectral type
    cat_tbl
        Pickles catalogue table: ``pyckles.SpectralLibrary("pickles").table``

    Returns
    -------
    new_spec : str
        Closest spectral type in the Pickles catalogue

    """
    if isinstance(spec, (list, tuple, np.ndarray)):
        return [nearest_spec_type(spt, cat_tbl) for spt in spec]
    else:
        if spec in cat_tbl["name"]:
            return spec

        elif spec[:2] in cat_tbl["luminosity"]:
            mask = np.where(cat_tbl["luminosity"] == spec[:2])[0]
            avail_evols = cat_tbl["evolution"][mask]
            roman = spec[4:] if spec[2] == "." else spec[2:]
            spec_evol = roman_to_arabic(roman)
            new_evol = avail_evols[np.argmin(np.abs(avail_evols - spec_evol))]
            new_spec = spec[:2] + arabic_to_roman(new_evol).upper()
            return new_spec

        elif spec[0] in "OBAFGKM":
            cat_lums = cat_tbl["luminosity"]
            mask = np.where([spec[0] == lum[0] for lum in cat_lums])[0]
            avail_lums = cat_tbl["luminosity"][mask]
            avail_class = np.array([spt[1] for spt in avail_lums]).astype(int)
            avail_class = np.unique(avail_class)
            spec_class = int(spec[1])
            new_class = avail_class[np.argmin(np.abs(avail_class - spec_class))]
            new_spec = spec[0] + str(new_class) + spec[2:]
            return nearest_spec_type(new_spec, cat_tbl)

        else:
            raise ValueError("{} was not even close to a Pickles entry"
                             "".format(spec))


def roman_to_arabic(x):
    dic = {"i": 1, "ii": 2, "iii": 3, "iv": 4, "v": 5, "vi": 6}
    return dic[x.lower()]


def arabic_to_roman(y):
    dic = {1: "i", 2: "ii", 3: "iii", 4: "iv", 5: "v", 6: "vi"}
    return dic[y]

