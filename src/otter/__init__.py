from __future__ import annotations
import logging
logger = logging.getLogger(__name__)

# get the version
from ._version import __version__


# explicitly set the package variable to ensure relative import work
__package__ = "otter"

# lazy load the SFD dust map, if needed
try:

    # this is a hacky fix for some header cleaning in dataverse
    # see https://github.com/gregreen/dustmaps/issues/72
    # TODO: This code block should be removed once that issue is figured out
    import requests.utils
    requests.utils.default_user_agent = lambda: "Mozilla/5.0"

    # now on with the normal dustmaps downloads :)
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
except ModuleNotFoundError:
    logger.warning(
        "Not loading dustmaps module! This means the photometry may or may not be \
        MW extinction corrected. We suggest installing dustmaps and/or checking the \
        photometry"
    )


# import important stuff
from .io.otter import Otter
from .io.transient import Transient
from . import util
from . import schema
from . import exceptions

# other, optional modules
try:
    from .io.host import Host
    from .io.data_finder import DataFinder
    from .plotter.otter_plotter import OtterPlotter
    from .plotter.plotter import plot_light_curve, plot_sed, quick_view, query_quick_view
except ModuleNotFoundError:
    logger.warning(
        "Not loading DataFinder, Host, and plotter modules for a minimal installation"
    )
