# -*- coding: utf-8 -*-
"""Store default values and keep ScopeSim imports in one place."""

from astar_utils import NestedMapping

import scopesim

Source = scopesim.source.source.Source
ter_curve_utils = scopesim.effects.ter_curves_utils
im_plane_utils = scopesim.optics.image_plane_utils
scopesim_utils = scopesim.utils
load_example_optical_train = scopesim.load_example_optical_train

__config__ = NestedMapping({
    "spectral": {
        "wave_min": 0.3,
        "wave_max": 3.0,
    },
    "random": {
        "seed": 9001,
    },
})
