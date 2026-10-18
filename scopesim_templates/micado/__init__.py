from astropy import units as u

from .darkness import darkness as _darkness
from .flatlamp import flatlamp as _flatlamp
from .pinhole_masks import pinhole_mask as _pinhole_mask
from .cluster import cluster as _cluster
from .viking_fields import viking_fields
from ..utils.general_utils import add_function_call_str

from . import spectral_calibrations

@add_function_call_str
def darkness():
    return _darkness()


@add_function_call_str
def empty_sky():
    # Include empty_sky for backwards compatibility reasons.
    # TODO: Remove empty_sky when it is not used in MicadoWISE anymore.
    return _darkness()


@add_function_call_str
def flatlamp():
    return _flatlamp()


@add_function_call_str
def pinhole_mask(x, y, waves, temperature=1500*u.K, sum_factor=1):
    return _pinhole_mask(x, y, waves, temperature=temperature,
                         sum_factor=sum_factor)


pinhole_mask.__doc__ = _pinhole_mask.__doc__


@add_function_call_str
def cluster():
    return _cluster()


__all__ = ["darkness", "flatlamp", "empty_sky", "cluster", "pinhole_mask",
           "viking_fields", "spectral_calibrations"]
