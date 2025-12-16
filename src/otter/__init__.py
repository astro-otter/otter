from __future__ import annotations

# get the version
from ._version import __version__


# explicitly set the package variable to ensure relative import work
__package__ = "otter"

# lazy load the SFD dust map, if needed
import os
import dustmaps
dustmaps_datapath = os.path.join(
    os.path.dirname(dustmaps.__file__),
    "data",
    "sfd",
    "SFD_dust_4096_sgp.fits"
)
if not os.path.exists(dustmaps_datapath):
    import dustmaps.sfd
    dustmaps.sfd.fetch()

# import important stuff
from .io.otter import Otter
from .io.transient import Transient
from .io.host import Host
from .io.data_finder import DataFinder
from .plotter.otter_plotter import OtterPlotter
from .plotter.plotter import plot_light_curve, plot_sed, quick_view, query_quick_view
from . import util
from . import schema
from . import exceptions
