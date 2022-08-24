from .darkness import darkness as _darkness
from .flatlamp import flatlamp as _flatlamp
from .cluster import cluster as _cluster
from .viking_fields import viking_fields
from ..utils.general_utils import add_function_call_str


@add_function_call_str
def darkness():
    _darkness()


@add_function_call_str
def empty_sky():
    # Include empty_sky for backwards compatibility reasons.
    # TODO: Remove empty_sky when it is not used in MicadoWISE anymore.
    return _darkness()


@add_function_call_str
def flatlamp():
    return _flatlamp()


@add_function_call_str
def cluster():
    return _cluster()


__all__ = [darkness, flatlamp, empty_sky, cluster]
