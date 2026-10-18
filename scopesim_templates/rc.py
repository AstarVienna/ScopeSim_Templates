# -*- coding: utf-8 -*-
"""Store default values and keep ScopeSim imports in one place."""

from astar_utils import NestedMapping

import scopesim

# Prefer ScopeSim's public API; the deep import is kept for older versions.
Source = getattr(scopesim, "Source", None) or scopesim.source.source.Source

# Not part of ScopeSim's public API. All ScopeSim internals used by this
# package are imported here, in one place, so a ScopeSim refactor only needs
# fixing in this module.
create_wcs_from_points = scopesim.optics.image_plane_utils.create_wcs_from_points
load_example_optical_train = scopesim.load_example_optical_train

__config__ = NestedMapping({
    "spectral": {
        "wave_min": 0.3,
        "wave_max": 3.0,
    },
    "random": {
        # Default seed for all stochastic templates. None means a fresh,
        # non-reproducible realisation on every call. Set an int here to
        # make every template reproducible without passing ``seed=``.
        "seed": None,
    },
})
