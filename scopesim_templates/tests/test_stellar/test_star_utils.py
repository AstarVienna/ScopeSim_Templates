import pytest

import pyckles

from scopesim_templates.stellar import stars_utils as stu


@pytest.fixture(scope="class")
def pickles():
    return pyckles.SpectralLibrary("pickles")


@pytest.mark.usefixtures("pickles")
class TestNearestSpecType:
    @pytest.mark.parametrize("spec_types", ("A0V", "G2V",
                                            "O8III", "G5III_R"))
    def test_returns_existing_spec_type(self, pickles, spec_types):
        spec = stu.nearest_spec_type(spec_types, pickles.table)
        assert spec == spec_types

    @pytest.mark.parametrize("input, output", (["K7I", "K7V"],
                                               ["G2II", "G2I"],
                                               ["G2III", "G2IV"]))
    def test_returns_nearest_evolutionary_type(self, pickles, input, output):
        spec = stu.nearest_spec_type(input, pickles.table)
        assert spec == output

    @pytest.mark.parametrize("input, output", (["O1I", "O5V"],
                                               ["G1II", "G0I"],
                                               ["F4V", "F5V"]))
    def test_returns_nearest_luminosity_class(self, pickles, input, output):
        spec = stu.nearest_spec_type(input, pickles.table)
        assert spec == output

    def test_returns_list_for_list_input(self, pickles):
        list_in = ["O1I", "G1II", "F4V"]
        list_out = ["O5V", "G0I", "F5V"]
        spec = stu.nearest_spec_type(list_in, pickles.table)
        assert all([so == lo for so, lo in zip(spec, list_out)])

    def test_throws_error_if_wierd_stuff_inputted(self, pickles):
        with pytest.raises(ValueError):
            stu.nearest_spec_type("bogus", pickles.table)
