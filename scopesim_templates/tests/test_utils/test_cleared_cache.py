from scopesim_templates.basic.stars import cluster
from scopesim import Source
from scopesim.effects.ter_curves_utils import download_svo_filter


def actually_works():
    my_cluster = cluster(mass=900, distance=50000, core_radius=1)
    assert isinstance(my_cluster, Source)
