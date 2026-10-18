# -*- coding: utf-8 -*-
"""Autouse fixture asserting that a test leaves the global NumPy RNG alone.

Import it into a test module with
``from ._rng_guard import global_rng_untouched  # noqa: F401``.
"""

import pytest
import numpy as np
from numpy import testing as npt


@pytest.fixture(autouse=True)
def global_rng_untouched():
    """Every test using this must leave the legacy global NumPy RNG alone."""
    state = np.random.get_state()
    yield
    after = np.random.get_state()
    assert state[0] == after[0]
    npt.assert_array_equal(state[1], after[1])
    assert state[2:] == after[2:]
