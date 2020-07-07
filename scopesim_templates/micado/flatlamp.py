"""Flat lamp for MICADO."""

import numpy as np
import scopesim


def flatlamp():
    """Create flat lamp source."""
    return scopesim.Source(
        x=[0], y=[0], ref=[0], spectra=np.array([1, 1]), lam=np.array([0.5, 2.5])
    )
