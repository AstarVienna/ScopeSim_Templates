"""A standard cluster for MICADO."""

from ..utils.general_utils import add_function_call_str
from ..stellar import clusters


@add_function_call_str
def cluster():
    """A basic cluster."""
    return clusters.cluster(mass=10000, distance=2000, core_radius=.1)
