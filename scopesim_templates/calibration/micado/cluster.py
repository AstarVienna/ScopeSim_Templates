"""A standard cluster for MICADO."""
import scopesim_templates.stellar.clusters


def cluster():
    """A basic cluster."""
    return scopesim_templates.stellar.clusters.cluster(mass=10000, distance=2000, core_radius=.1)
