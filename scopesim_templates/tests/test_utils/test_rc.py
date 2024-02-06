from astar_utils import NestedMapping

from scopesim_templates import rc


class TestRC:
    def test_source_object_is_imported(self):
        assert hasattr(rc, "Source")

    def test_default_yaml_is_read_in(self):
        isinstance(rc.__config__, NestedMapping)
