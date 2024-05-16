"""
Some utilities to create common plots for transients that use the OtterPlotter
"""

from __future__ import annotations
from .otter_plotter import OtterPlotter


def plot_light_curve(
    date: float,
    flux: float,
    date_err: float = None,
    flux_err: float = None,
    fig=None,
    ax=None,
    backend: str = "matplotlib",
    xlabel: str = "Date",
    ylabel: str = "Flux",
    **kwargs,
):
    """
    Plot the light curve for the input data

    Args:
        date (float): MJD dates
        flux (float): Flux
        date_err (float): optional error on the MJD dates
        flux_err (float): optional error on the flux
        fig (float): matplotlib fig object, optional. Will be created if not provided.
        ax (float): matplitlib axis object, optional. Will be created if not provided.
        backend (str): backend for plotting. options: "matplotlib" (default) or "plotly"
        xlabel (str): x-axis label
        ylabel (str): y-axis label
        **kwargs: keyword arguments to pass to either plotly.graph_objects.add_scatter
                  or matplotlib.pyplot.errorbar

    Returns:
       Either a matplotlib axis or plotly figure
    """

    plt = OtterPlotter(backend)
    fig = plt.plot(date, flux, date_err, flux_err, ax=ax, **kwargs)

    if backend == "matplotlib":
        fig.set_ylabel(ylabel)
        fig.set_xlabel(xlabel)

    elif backend == "plotly":
        fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel, **kwargs)

    return fig


def plot_sed(
    wave_or_freq: float,
    flux: float,
    wave_or_freq_err: float = None,
    flux_err: float = None,
    fig=None,
    ax=None,
    backend: str = "matplotlib",
    xlabel: str = "Frequency or Wavelength",
    ylabel: str = "Flux",
    **kwargs,
):
    """
    Plot the SED for the input data

    Args:
        wave_or_freq (float): wave or frequency array
        flux (float): Flux
        wave_or_freq_err (float): optional error on the MJD dates
        flux_err (float): optional error on the flux
        fig (float): matplotlib fig object, optional. Will be created if not provided.
        ax (float): matplitlib axis object, optional. Will be created if not provided.
        backend (str): backend for plotting. Options: "matplotlib" (default) or "plotly"
        xlabel (str): x-axis label
        ylabel (str): y-axis label
        **kwargs: keyword arguments to pass to either plotly.graph_objects.add_scatter
                  or matplotlib.pyplot.errorbar

    Returns:
       Either a matplotlib axis or plotly figure
    """

    plt = OtterPlotter(backend)
    fig = plt.plot(wave_or_freq, flux, wave_or_freq_err, flux_err, ax=ax, **kwargs)

    if backend == "matplotlib":
        fig.set_ylabel(ylabel)
        fig.set_xlabel(xlabel)

    elif backend == "plotly":
        fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel, **kwargs)

    return fig
