from pathlib import Path
import sys
import yaml

from astar_utils import NestedMapping

# TODO: what's this?
if Path("F:/Work/ScopeSim/").exists():
    sys.path.append("F:/Work/ScopeSim/")

import scopesim

Source = scopesim.source.source.Source
ter_curve_utils = scopesim.effects.ter_curves_utils
im_plane_utils = scopesim.optics.image_plane_utils
scopesim_utils = scopesim.utils
load_example_optical_train = scopesim.load_example_optical_train

PKG_DIR = Path(__file__).parent
with (PKG_DIR / "defaults.yaml").open(encoding="utf-8") as file:
    __config__ = NestedMapping(yaml.full_load(file))
