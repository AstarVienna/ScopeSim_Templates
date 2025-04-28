# -*- coding: utf-8 -*-
"""."""

import numpy as np
from astropy import units as u

from spextra import Spextrum

from ..rc import Source
from . import source_from_array


def ellipse_2d(
        a_outer: float | int,
        a_inner: float | int,
        axis_ratio: float,
        angle: u.Quantity[u.deg] | float,
        width: int,
        height: int
) -> np.ndarray:
    """
    Draw an ellipse with a flux profile into a 2D array.

    Parameters
    ----------
    a_outer : float | int
        Donut ellipse outer semi-major axis.
    a_inner : float | int
        Donut ellipse inner semi-major axis.
    axis_ratio : float
        Ratio of minor to major axis 0.0-1.0.
    angle : u.Quantity[u.deg] | float
        Orientation of the ellipse's major axis. If float assumes deg.
        Zero aligns major axis with the x-axis.
    width : int
        Image width in pixel.
    height : int
        Image heigth in pixel.

    Returns
    -------
    image : np.ndarray
        Output image.

    """
    cen = (np.array([width, height]) - 1) / 2
    angle = (angle << u.deg << u.rad).value

    x_grid, y_grid = np.meshgrid(np.arange(width) - cen[0],
                                 np.arange(height) - cen[1])
    # Rotate the coordinate grids
    xp_grid = +x_grid * np.cos(angle) + y_grid * np.sin(angle)
    yp_grid = -x_grid * np.sin(angle) + y_grid * np.cos(angle)

    # Elliptical radius squared
    rell2_grid = xp_grid**2 + yp_grid**2 / axis_ratio**2

    # Define the donut as a mask
    donutmask = ((rell2_grid <= a_outer**2) &
                 (rell2_grid >= a_inner**2))

    # Brightness profile - this is very arbitrary
    img = (50 + 300 * np.exp(-rell2_grid / (2 * (0.2 * width)**2))) * donutmask
    return img


def donut(
        outer=1*u.arcsec,
        fraction=.4,
        inclination=30*u.deg,
        angle=30*u.deg,
        amplitude=16*u.ABmag,
        filter_curve="R",
        sed=None,
        pixel_scale=0.00057*u.arcsec,  # where does the number come from?
        width=None,
        height=None,
) -> Source:
    """Define an disk with a central hole (a "donut") seen at an inclination

    Parameters
    ----------
    outer : u.Quantity[u.arcsec] | float
        Donut outer radius. If float assumes arcsec.
    fraction : float
        Filling fraction of the disk. The inner radius is (1 - fraction) * outer.
    inclination : u.Quantity[u.deg] | float
        Inclination angle of the disk (0: face on; 90 deg: edge-on)
        If float assumes deg.
    angle : u.Quantity[u.deg] | float
        Orientation of the inclined donut's major axis. If float assumes deg.
        Zero aligns major axis with the x-axis.
    """
    if sed is None:
        sed = Spextrum.black_body_spectrum(2000*u.K, amplitude, filter_curve)
    outer = (outer << u.arcsec).value
    pixel_scale = (pixel_scale << u.arcsec).value
    outer = outer / pixel_scale
    inner = (1 - fraction) * outer
    inclination = (inclination << u.deg << u.rad).value
    axis_ratio = np.cos(inclination)

    # Size of the output image is computed from the extent of the donut
    if width is None:
        width = (2.1 * outer).astype(int)
    if height is None:
        height = (2.1 * outer).astype(int)

    img = ellipse_2d(outer, inner, axis_ratio, angle, width, height)
    src = source_from_array(
        arr=img, sed=sed,
        pixel_scale=pixel_scale,
        amplitude=amplitude,
        filter_curve=filter_curve,
    )
    return src
