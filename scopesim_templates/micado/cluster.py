"""A standard cluster for MICADO."""

from scopesim_templates import basic


def cluster():
    """A basic cluster."""
    return basic.stars.cluster(mass=10000, distance=2000, core_radius=.1)
