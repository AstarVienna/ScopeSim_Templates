import numpy as np

from ..rc import Source, __config__


def empty_sky():
    """
    Returns an empty source so that instrumental fluxes can be simulated

    Returns
    -------
    sky : Source

    """
    sky = Source(lam=np.array([__config__["!spectral.wave_min"],
                               __config__["!spectral.wave_max"]]),
                 spectra=np.array([0, 0]), x=[0], y=[0], ref=[0], weight=[0])
    return sky
