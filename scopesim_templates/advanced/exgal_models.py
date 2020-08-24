"""
This file should contain analytical extragalactic models
"""


import numpy as np

from astropy import units as u
from astropy.units import Quantity, UnitsError
from astropy.modeling.core import (Fittable1DModel, Fittable2DModel,
                                   ModelDefinitionError)

from astropy.modeling.parameters import Parameter, InputParameterError
from astropy.modeling.models import Sersic2D


# Definitions

TWOPI = 2 * np.pi
FLOAT_EPSILON = float(np.finfo(np.float32).tiny)
GAUSSIAN_SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


class VelField(Fittable2DModel):
    r"""
        Two dimensional Velocity Field following arctan approximation

        Parameters
        ----------

        ellip : float, u.Quantity
            Ellipticity on the sky

        theta : float, u.Quantity
            Position angle of the major axis wrt to north (=up) measured counterclockwise,
        vmax : float, u.Quantity
            Constant rotation velocity for R>>rd,

        r_eff : float
            scale length of galaxy (assumed to be turnover radius)

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

    @staticmethod
    def evaluate(x, y, vmax, r_eff,  ellip, theta, x_0, y_0, q):
        """
        Two dimensional velocity field, arctan approximation

        """
        if isinstance(theta, u.Quantity) is False:
            theta = theta * u.deg

        r_d = r_eff  # For now,  for n=1  r_eff = 1.678 * r_d
        theta = (-theta).to(u.rad)
        # get inclination from ellipticity
        incl = np.arccos(np.sqrt(((1 - ellip) ** 2 - q ** 2) / (1 - q ** 2)))

        r = ((x - x_0) ** 2 + (y - y_0) ** 2) ** 0.5

        #   azimuthal angle in the plane of the galaxy = cos(theta) = cost
        cost = (-(x - x_0) * np.sin(theta) + (y - y_0) * np.cos(theta)) / (r + 0.00001)
        vrot = vmax*2 / np.pi*np.arctan(r/r_d)         #arctan model

        return vrot * np.sin(incl) * cost

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        else:
            return {'x': self.x_0.unit,
                    'y': self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit['x'] != inputs_unit['y']:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit['x'],
                'y_0': inputs_unit['x'],
                'r_eff': inputs_unit['x'],
                'phi': u.deg,
                'vrot': outputs_unit['z']}



class DispersionField(Fittable2DModel):

    r"""
        Two dimensional Velocity Dispersion
        This follows a broken power law as in Veale et al. 2017

        .. math:: \sigma(R)=\sigma_0 2^{\gamma_1-\gamma_2} (\frac{R}{R_b})^{\gamma_1}(1+\frac{R}{R_b}^{\gamma_2-\gamma_1}

        where
        :math:`\sigma_0` is the velocity dispersion
        :math:`\gamma_1`  and math:`\gamma_2` are the inner and outer power slopes. Set to -0.04 and -0.42 as
        default values
        :math:`R_b` is the break radius, set to :math:`R_b = R_{eff}` for simplicity


        Parameters
        ----------

        incl : float, u.Quantity
            Inclination inclination between the normal to the galaxy plane and the line-of-sight,

        phi : float, u.Quantity
            Position angle of the major axis wrt to north (=up) measured counterclockwise,
        sigma : float, u.Quantity
            velocity dispersion

        r_eff : float
            scale length of galaxy (assumed to be turnover radius) it will be used as sigma

        x0 : float, optional
            x position of the center.
        y0 : float, optional
            y position of the center.

        e_in : float, optional
            inner power law slope

        e_out : float, optional
            outer power law slope

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

    @staticmethod
    def evaluate(x, y, ellip, theta, sigma, r_eff, x_0, y_0, e_in, e_out):
        """
        evaluation

        """

        if isinstance(theta, u.Quantity) is False:
            theta = theta * u.deg

        theta = theta.to(u.rad)

        a, b = r_eff, (1 - ellip) * r_eff
        cos_theta, sin_theta = np.cos(theta), np.sin(theta)
        x_maj = (x - x_0) * sin_theta + (y - y_0) * cos_theta + 0.1
        x_min = -(x - x_0) * cos_theta + (y - y_0) * sin_theta + 0.1  # to avoid inf values
        z = np.sqrt((x_maj / a) ** 2 + (x_min / b) ** 2)
        result = sigma * 2**(e_in - e_out) * z**e_in * (1 + z)**(e_out - e_in)

        return result



    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        else:
            return {'x': self.x_0.unit,
                    'y': self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit['x'] != inputs_unit['y']:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return {'x_0': inputs_unit['x'],
                'y_0': inputs_unit['x'],
                'r_eff': inputs_unit['x'],
                'phi': u.deg,
                'sigma': outputs_unit['z']}


class GalaxyBase:
    """
    Class to define a galaxy. It takes a set of parameters and creates the
    different moments for

    """


    def __init__(self, x, y, x_0, y_0,
                 amplitude, r_eff, ellip, theta, n=4,
                 vmax=0, sigma=0, q=0.2,
                 e_in=-0.04, e_out=-0.42):

        self.x = x
        self.y = y
        self.amplitude = amplitude
        self.r_eff = r_eff
        self.x_0 = x_0
        self.y_0 = y_0
        self.ellip = ellip
        self.theta = theta
        self.n = n
        self.vmax = vmax
        self.sigma = sigma
        self.q = q
        self.e_in = e_in
        self.e_out = e_out

    @property
    def intensity(self):
        """
        2D light distribution following a Sersic2D profile

        Returns
        -------
        numpy array
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
    def velocity(self):
        """
        Velocity field according to the supplied parameters

        Returns
        -------
        numpy array
        """
        if self.vmax > 0:
            mod = VelField(x_0=self.x_0,
                           y_0=self.y_0,
                           r_eff=self.r_eff,
                           ellip=self.ellip,
                           theta=self.theta,
                           vmax=self.vmax,
                           q=self.q)
            result = mod(self.x, self.y)
        else:
            result = np.ones(shape=self.x.shape)

        return result

    @property
    def dispersion(self):
        """
        Velocity dispersion map according to the supplied parameters

        Returns
        -------
        numpy array
        """
        if self.sigma > 0:
            mod = DispersionField(x_0=self.x_0,
                                  y_0=self.y_0,
                                  r_eff=self.r_eff,
                                  ellip=self.ellip,
                                  theta=self.theta,
                                  sigma=self.sigma,
                                  e_in=self.e_in,
                                  e_out=self.e_out)
            result = mod(self.x, self.y)
        else:
            result = np.ones(shape=self.x.shape)

        return result

    def regrid(self, ngrid=10):
        """
        Regrid the smooth velocity field to regions with similar
        velocity and velocity dispersion

        Parameters
        ----------
        ngrid: integer

        Returns
        -------
        A numpy array with sectors numbered
        """
        velfield = self.velocity.value
        dispfield = self.dispersion.value

        vel_grid = np.round((ngrid // 2) * velfield / np.max(velfield)) * np.max(velfield)
        sigma_grid = np.round((ngrid // 2 + 2) * dispfield / np.max(dispfield)) * np.max(dispfield)
        total_field = vel_grid + sigma_grid
        uniques = np.unique(total_field)

        for i, v in enumerate(uniques):
            total_field[total_field == v] = i+1

        return total_field

    def get_masks(self, ngrid=10):
        """
        Returns the regrided regions as a list of numpy masks

        Parameters
        ----------
        ngrid: integer, griding factor

        Returns
        -------
        list of masks
        """

        grid = self.regrid(ngrid=ngrid)
        uniques = np.unique(grid)
        masklist = []
        for value in uniques:
            mask = np.ma.masked_where(grid == value, grid, copy=True).mask.astype(int)
            masklist.append(mask)

        return masklist


