import warnings

from . import calibration
from . import extragalactic
from . import misc
from . import stellar
from . import utils
from . import metis

# This warning is emitted when just doing "import scopesim_templates", which
# should be a normal non-warning thing to do.
# TODO: Find a way to emit this warning only when the functions below are
#       accessed directly.
# warnings.warn("In a future version top level function calls will be removed. "
#               "Always use this syntax: from module.submodule import function",
#               DeprecationWarning, stacklevel=2)
from .misc.misc import source_from_image, source_from_file, source_from_cube
from .stellar.stars import star, stars, star_grid, star_field
from .stellar.clusters import cluster
from .extragalactic.galaxies import galaxy, galaxy3d, elliptical
from .calibration.calibration import empty_sky
from .micado import darkness, flatlamp
from . import micado
