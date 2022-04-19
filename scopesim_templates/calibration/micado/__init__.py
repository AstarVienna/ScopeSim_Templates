from .darkness import darkness
from .flatlamp import flatlamp
from .cluster import cluster
from .viking_fields import viking_fields
# Include empty_sky for backwards compatibility reasons.
# TODO: Remove empty_sky when it is not used in MicadoWISE anymore.
empty_sky = darkness

__all__ = [darkness, flatlamp, empty_sky, cluster]
