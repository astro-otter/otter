from __future__ import annotations

# get the version
from ._version import __version__


# explicitly set the package variable to ensure relative import work
__package__ = "otter"

# import important stuff
from .io.otter import Otter
from .io.transient import Transient
from .plotter.otter_plotter import OtterPlotter
from .plotter.plotter import plot_light_curve, plot_sed
