"""
OtterPlotter Class which handles the backend of the plotting.

Currently supported backends are:
- matplotlib
- plotly
"""

from __future__ import annotations
import importlib


class OtterPlotter:
    """
    Handles the backend for the "plotter" module

    Args:
        backend (string): a string of the module name to import and use
                          as the backend. Currently supported are "matplotlib",
                          "matplotlib.pyplot", "plotly", and "plotly.graph_objects"
    """

    def __init__(self, backend):
        if backend == "matplotlib.pyplot":
            self.backend = backend
        elif backend == "pyplot.graph_objects":
            self.backend = backend
        elif "plotly" in backend and "graph_objects" not in backend:
            self.backend = "plotly.graph_objects"
        elif "matplotlib" in backend and "pyplot" not in backend:
            self.backend = "matplotlib.pyplot"
        else:
            raise ValueError("Not a valid backend string!")

        self.plotter = importlib.import_module(self.backend)

        if self.backend == "matplotlib.pyplot":
            self.plot = self._plot_matplotlib
        elif self.backend == "plotly.graph_objects":
            self.plot = self._plot_plotly
        else:
            raise ValueError("Unknown plotting backend!")

    def _plot_matplotlib(self, x, y, xerr=None, yerr=None, ax=None, **kwargs):
        """
        General plots using matplotlib, is called by _matplotlib_light_curve and
        _matplotlib_sed
        """

        if ax is None:
            _, ax = self.plotter.subplots()

        ax.errorbar(x, y, xerr=xerr, yerr=yerr, **kwargs)
        return ax

    def _plot_plotly(self, x, y, xerr=None, yerr=None, ax=None, *args, **kwargs):
        """
        General plotting method using plotly, is called by _plotly_light_curve and
        _plotly_sed
        """

        if ax is None:
            go = self.plotter.Figure()
        else:
            go = ax

        fig = go.add_scatter(
            x=x, y=y, error_x=dict(array=xerr), error_y=dict(array=yerr), **kwargs
        )

        return fig
