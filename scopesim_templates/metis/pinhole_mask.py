# -*- coding: utf-8 -*-
"""Pinhole mask for METIS."""

from itertools import product

import numpy as np
from astropy import units as u

from ..utils.general_utils import add_function_call_str
from ..micado.pinhole_masks import pinhole_mask as micado_mask


@add_function_call_str
def pinhole_mask(
        nx=5,
        ny=28,
        dx=0.45,
        dy=0.27945,
        waves=np.arange(3, 6, .001)*u.um,
        temperature = 1500*u.K,
        **kwargs
):
    """
    Create pinhole mask source for METIS, simulated as a grid of point sources.

    This currently uses the (flexible)  MICADO pinhole mask internally.

    Parameters
    ----------
    nx : int, optional
        Number of points along x axis. The default is 5.
    ny : int, optional
        Number of points along y axis. The default is 28.
    dx : float, optional
        (Half) extension in x direction. The default is 0.45.
    dy : float, optional
        (Half) extension in y direction. The default is 0.27945.
    waves : array-like, optional
        Wavelength axis. The default is np.arange(3, 6, .001)*u.um.
    **kwargs
        Any other kwargs passed to micado.pinhole_mask.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    x, y = zip(*product(np.linspace(-dx, dx, nx), np.linspace(-dy, dy, ny)))

    return micado_mask(x, y, waves, temperature, **kwargs)

