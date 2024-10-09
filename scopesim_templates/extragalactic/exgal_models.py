"""Module should contain analytical extragalactic models."""

from dataclasses import dataclass
from collections.abc import Generator

import numpy as np

from astropy import units as u
from astropy.units import UnitsError
from astropy.modeling.core import Fittable2DModel
from astropy.modeling.parameters import Parameter
from astropy.modeling.models import Sersic2D


# Definitions

TWOPI = 2 * np.pi
FLOAT_EPSILON = float(np.finfo(np.float32).tiny)
GAUSSIAN_SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


class VelField(Fittable2DModel):
    r"""
    Two dimensional Velocity Field following arctan approximation.

    Parameters
    ----------
    ellip : float, u.Quantity
        Ellipticity on the sky
    theta : float, u.Quantity
        Position angle of the major axis. The rotation angle increases
        counterclockwise from the positive x axis.
    vmax : float, u.Quantity
        Constant rotation velocity for R>>rd.
    r_eff : float
        scale length of galaxy (assumed to be turnover radius).
    x0 : float, optional
        x position of the center.
    y0 : float, optional
        y position of the center.
    q : float, optional
        Disk thickness
    """

    vmax = Parameter(default=100)
    r_eff = Parameter(default=1)

    ellip = Parameter(default=0)    # maximum ellipticity 1 - q,  make tests
    theta = Parameter(default=0)

    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)

    q = Parameter(default=0.2)

    @classmethod
    def from_galaxy_base(cls, galaxy_base):
        """Create instance from ``GalaxyBase`` instance."""
        return cls(x_0=galaxy_base.x_0,
                   y_0=galaxy_base.y_0,
                   r_eff=galaxy_base.r_eff,
                   ellip=galaxy_base.ellip,
                   theta=galaxy_base.theta,
                   vmax=galaxy_base.vmax,
                   q=galaxy_base.q)

    @staticmethod
    def evaluate(x, y, vmax, r_eff,  ellip, theta, x_0, y_0, q):
        """Two dimensional velocity field, arctan approximation."""
        theta <<= u.deg
        theta = (-theta).to(u.rad)

        r_d = r_eff  # For now,  for n=1  r_eff = 1.678 * r_d
        # get inclination from ellipticity
        incl = np.arccos(np.sqrt(((1 - ellip) ** 2 - q ** 2) / (1 - q ** 2)))

        r = ((x - x_0) ** 2 + (y - y_0) ** 2) ** 0.5

        #   azimuthal angle in the plane of the galaxy = cos(theta) = cost
        cost = (-(x - x_0) * np.cos(theta) + (y - y_0) *
                np.sin(theta)) / (r + 0.00001)
        vrot = vmax*2 / np.pi*np.arctan(r/r_d)         # arctan model

        return vrot * np.sin(incl) * cost

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {"x": self.x_0.unit,
                "y": self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit["x"] != inputs_unit["y"]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {"x_0": inputs_unit["x"],
                "y_0": inputs_unit["x"],
                "r_eff": inputs_unit["x"],
                "phi": u.deg,
                "vrot": outputs_unit["z"]}


class DispersionField(Fittable2DModel):
    r"""
    Two dimensional Velocity Dispersion.

    This follows a broken power law as in Veale et al. 2017

    .. math:: \sigma(R)=\sigma_0 2^{\gamma_1-\gamma_2} \left(\frac{R}{R_b}\right)^{\gamma_1}\left(1+\frac{R}{R_b}\right)^{\gamma_2-\gamma_1}

    where
    :math:`\sigma_0` is the velocity dispersion
    :math:`\gamma_1` and :math:`\gamma_2` are the inner and outer power slopes.
    Set to -0.04 and -0.42 as default values
    :math:`R_b` is the break radius, set to :math:`R_b = R_{eff}` for
    simplicity


    Parameters
    ----------
    incl : float, u.Quantity
        Inclination inclination between the normal to the galaxy plane and the
        line-of-sight.
    phi : float, u.Quantity
        Position angle of the major axis wrt to north (=up) measured
        counterclockwise.
    sigma : float, u.Quantity
        Velocity dispersion
    r_eff : float
        Scale length of galaxy (assumed to be turnover radius) it will be used
        as sigma.
    x0 : float, optional
        x position of the center.
    y0 : float, optional
        y position of the center.
    e_in : float, optional
        Inner power law slope.
    e_out : float, optional
        Outer power law slope.
    theta : float, u.Quantity
        Position angle of the major axis. The rotation angle increases
        counterclockwise from the positive x axis.
    """

    # incl = Parameter(default=45)
    ellip = Parameter(default=0)
    theta = Parameter(default=0)
    sigma = Parameter(default=100)
    r_eff = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    e_in = Parameter(default=-0.0)
    e_out = Parameter(default=-0.5)

    @classmethod
    def from_galaxy_base(cls, galaxy_base):
        """Create instance from ``GalaxyBase`` instance."""
        return cls(x_0=galaxy_base.x_0,
                   y_0=galaxy_base.y_0,
                   r_eff=galaxy_base.r_eff,
                   ellip=galaxy_base.ellip,
                   theta=galaxy_base.theta,
                   sigma=galaxy_base.sigma,
                   e_in=galaxy_base.e_in,
                   e_out=galaxy_base.e_out)

    @staticmethod
    def evaluate(x, y, ellip, theta, sigma, r_eff, x_0, y_0, e_in, e_out):
        """Perform evaluation."""
        theta <<= u.deg
        theta = theta.to(u.rad)

        a, b = r_eff, (1 - ellip) * r_eff
        sin_theta, cos_theta = np.sin(theta), np.cos(theta)
        x_maj = (x - x_0) * cos_theta + (y - y_0) * sin_theta + 0.1
        x_min = -(x - x_0) * sin_theta + (y - y_0) * cos_theta + 0.1  # to avoid inf values
        z = np.sqrt((x_maj / a) ** 2 + (x_min / b) ** 2)
        result = sigma * 2**(e_in - e_out) * z**e_in * (1 + z)**(e_out - e_in)

        return result

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        return {"x": self.x_0.unit,
                "y": self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit["x"] != inputs_unit["y"]:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {"x_0": inputs_unit["x"],
                "y_0": inputs_unit["x"],
                "r_eff": inputs_unit["x"],
                "phi": u.deg,
                "sigma": outputs_unit["z"]}


@dataclass
class GalaxyBase:
    """
    Class to define a galaxy.

    It takes a set of parameters and creates the different moments for
    """

    x: np.ndarray
    y: np.ndarray
    x_0: int
    y_0: int
    amplitude: float
    r_eff: float
    ellip: float
    theta: float
    n: int = 4
    vmax: float = 0.0
    sigma: float = 0.0
    q: float = 0.2
    e_in: float = 0.04
    e_out: float = -0.42

    @property
    def intensity(self) -> np.ndarray:
        """
        2D light distribution following a Sersic2D profile.

        Returns
        -------
        result : numpy.ndarray
            Light distribution map.

        """
        mod = Sersic2D(x_0=self.x_0,
                       y_0=self.y_0,
                       amplitude=self.amplitude,
                       r_eff=self.r_eff,
                       n=self.n,
                       ellip=self.ellip,
                       theta=self.theta)
        return mod(self.x, self.y)

    @property
    def velocity(self) -> np.ndarray:
        """
        Velocity field according to the supplied parameters.

        Returns
        -------
        result : numpy.ndarray
            Velocity map.

        """
        if self.vmax <= 0.0:
            return np.ones_like(self.x)
        mod = VelField.from_galaxy_base(self)
        return mod(self.x, self.y)

    @property
    def dispersion(self) -> np.ndarray:
        """
        Velocity dispersion map according to the supplied parameters.

        Returns
        -------
        result : numpy.ndarray
            Velocity dispersion map.

        """
        if self.sigma <= 0.0:
            return np.ones_like(self.x)
        mod = DispersionField.from_galaxy_base(self)
        return mod(self.x, self.y)

    def regrid(self, ngrid: int = 10) -> np.ndarray:
        """
        Regrid the smooth velocity field.

        Regrid to regions with similar velocity and velocity dispersion.

        Parameters
        ----------
        ngrid : int, optional
            Griding factor. The default is 10.

        Returns
        -------
        sectors : numpy.ndarray
            A numpy array with sectors numbered.
        uniques : set
            Set of sector IDs.

        """
        velfield = self.velocity.value
        dispfield = self.dispersion.value

        vel_grid = np.round((ngrid // 2) * velfield / velfield.max())
        # Get order of magnitude for offset
        offset = 10**(int(np.log10(vel_grid.max())) + 2)
        sigma_grid = np.round((ngrid // 2 + 2) * dispfield /
                              dispfield.max()) * offset
        total_field = vel_grid + sigma_grid
        _, sectors = np.unique(total_field, return_inverse=True)
        uniques = set(sectors)
        # logger.debug("%d sectors", len(uniques))

        return sectors.reshape(total_field.shape), uniques

    def get_masks(self, ngrid: int = 10) -> Generator[np.ndarray]:
        """
        Return a generator of numpy masks from the regrided regions.

        Parameters
        ----------
        ngrid : TYPE, optional
            Griding factor, passed to ``GalaxyBase.regrid()``.
            The default is 10.

        Yields
        ------
        mask : numpy.ndarray
            Boolean array constructed from regrided regions.

        """
        grid, uniques = self.regrid(ngrid=ngrid)
        for value in uniques:
            yield np.equal(grid, value)
