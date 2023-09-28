from scopesim_templates import cluster
from scopesim_templates.rc import Source


def test_actually_works():
    my_cluster = cluster(mass=900, distance=50000, core_radius=1)
    assert isinstance(my_cluster, Source)
