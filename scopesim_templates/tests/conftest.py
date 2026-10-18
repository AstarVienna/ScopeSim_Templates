"""Shared pytest configuration for the scopesim_templates test suite."""

import numpy as np
import pytest


def _global_rng_state():
    state = np.random.get_state(legacy=False)
    return state["state"]["key"].copy(), state["state"]["pos"]


@pytest.fixture(autouse=True)
def global_rng_is_untouched(request):
    """
    Fail any test after which NumPy's *global* random state has changed.

    Templates must draw random numbers from their own ``np.random.Generator``
    (see ``utils.general_utils.get_rng``) and never reseed or consume the
    global stream, which belongs to the user. Tests that run a full ScopeSim
    observation may legitimately consume it (e.g. detector noise); mark those
    with ``@pytest.mark.allow_global_rng``.
    """
    if request.node.get_closest_marker("allow_global_rng"):
        yield
        return

    key_before, pos_before = _global_rng_state()
    yield
    key_after, pos_after = _global_rng_state()
    assert pos_before == pos_after and np.array_equal(key_before, key_after), (
        "NumPy's global random state was modified during this test")
