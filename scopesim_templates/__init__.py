import warnings

from . import calibration
from . import extragalactic
from . import misc
from . import stellar
from . import utils

warnings.warn("In a future version top level function calls will be removed. "
              "Always use this syntax: from module.submodule import function",
              DeprecationWarning)
from .misc.misc import source_from_image, source_from_file, source_from_cube
from .stellar.stars import star, stars, star_grid, star_field
from .stellar.clusters import cluster
from .extragalactic.galaxies import galaxy, galaxy3d, elliptical
from .calibration.calibration import empty_sky
from .calibration.micado import darkness, flatlamp
from .calibration import micado