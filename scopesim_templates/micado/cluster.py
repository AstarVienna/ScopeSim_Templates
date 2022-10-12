"""A standard cluster for MICADO."""
from scopesim_templates.utils.general_utils import function_call_str, add_function_call_str

import scopesim_templates.stellar.clusters


@add_function_call_str
def cluster():
    """A basic cluster."""
    return scopesim_templates.stellar.clusters.cluster(mass=10000, distance=2000, core_radius=.1)
