import numpy as np
from astropy.modeling.models import Sersic2D


# TODO: could this be done with the existing implementation of GalaxyBase?
def sersic_profile(r_eff=100, n=4, ellipticity=0.5, angle=30,
                   normalization="total",
                   width=1024, height=1024, x_offset=0, y_offset=0,
                   oversample=1):
    """
    Return a 2D array with a normalised Sersic profile.

    Parameters
    ----------
    r_eff : float
        [pixel] Effective radius

    n : float
        Power law index.
        - n=1 for exponential (spiral),
        - n=4 for de Vaucouleurs (elliptical)

    ellipticity : float
        Ellipticity is defined as (a - b)/a. Default = 0.5

    angle : float
        [deg] Default = 30. Rotation anti-clockwise from the x-axis

    normalization : str, optional
        ["half-light", "centre", "total"] Where the profile equals unity
        If normalization equals:
        - "half-light" : the pixels at the half-light radius are set to 1
        - "centre" : the maximum values are set to 1
        - "total" : the image sums to 1

    width, height : int
        [pixel] Dimensions of the image

    x_offset, y_offset : float
        [pixel] The distance between the centre of the profile and the centre
        of the image

    oversample : int
        Factor of oversampling, default factor = 1. If > 1, the model is
        discretized by taking the average of an oversampled grid.


    Returns
    -------
    img : 2D array


    Notes
    -----
    Most units are in [pixel] in this function. This differs from
    :func:`.galaxy` where parameter units are in [arcsec] or [pc]

    """
    # Silently cast to integer
    os_factor = int(oversample)

    if os_factor <= 0:
        raise ValueError("Oversampling factor must be >=1.")

    width_os = os_factor * width
    height_os = os_factor * height
    x, y = np.meshgrid(np.arange(width_os), np.arange(height_os))

    dx = 0.5 * width_os + x_offset * os_factor
    dy = 0.5 * height_os + y_offset * os_factor

    r_eff_os = r_eff * os_factor

    mod = Sersic2D(amplitude=1, r_eff=r_eff_os, n=n, x_0=dx, y_0=dy,
                   ellip=ellipticity, theta=np.deg2rad(angle))
    img_os = mod(x, y)

    # Rebin os_factord image
    img = _rebin(img_os, os_factor)

    thresh = np.max([img[0, :].max(), img[-1, :].max(),
                     img[:, 0].max(), img[:, -1].max()])
    img[img < thresh] = 0

    if "cen" in normalization.lower():
        img /= np.max(img)
    elif "tot" in normalization.lower():
        img /= np.sum(img)

    return img


def _rebin(img, bpix):
    """Rebin image img by block averaging bpix x bpix pixels."""
    xedge = np.shape(img)[0] % bpix
    yedge = np.shape(img)[1] % bpix
    img_block = img[xedge:, yedge:]

    binim = np.reshape(img_block,
                       (int(img_block.shape[0]/bpix), bpix,
                        int(img_block.shape[1]/bpix), bpix))
    binim = np.mean(binim, axis=3)
    binim = np.mean(binim, axis=1)
    return binim
