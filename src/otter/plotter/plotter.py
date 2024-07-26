"""
Some utilities to create common plots for transients that use the OtterPlotter
"""

from __future__ import annotations
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .otter_plotter import OtterPlotter
from ..exceptions import FailedQueryError
from ..io.otter import Transient, Otter


def query_quick_view(
    db: Otter = None,
    otter_path: str = None,
    ptype: str = "both",
    sed_dim: str = "freq",
    dt_over_t: float = 0,
    plotting_kwargs: dict = {},
    phot_cleaning_kwargs: dict = {},
    result_length_tol=10,
    **kwargs,
) -> list[plt.Figure]:
    """
    Queries otter and then plots all of the transients. It will either query the Otter
    object provided in db or construct an Otter object from otter_path.

    Args:
        db (otter.Otter) : The otter object to query
        otter_path (str) : The path to construct an otter path from
        ptype (str) : The plot type to generate. Valid options are
                      - both -> Plot both light curve and sed (default)
                      - sed -> Plot just the sed
                      - lc -> Plot just the light curve
        sed_dim (str) : The x dimension to plot in the SED. Options are "freq" or
                        "wave". Default is "freq".
        time_tol (float) : The tolerance to split the days by. Default is 1 day. must be
                           in units of days.
        plotting_kwargs (dict) : dictionary of key word arguments to pass to
                                 otter.plotter.plot_light_curve or
                                 otter.plotter.plot_sed.
        phot_cleaning_kwargs (dict) : Keyword arguments passed to
                                      otter.Transient.clean_photometry
        result_length_tol (int) : If the query result is longer than this it will throw
                                  an erorr to prevent 100s of plots from spitting out
                                  (and likely crashing your computer). Default is 10 but
                                  can be adjusted.
        **kwargs : Arguments to pass to otter.Otter.query

    Returns:
        A list of matplotlib pyplot Figure objects that we plotted

    """
    if db is None:
        if otter_path is not None:
            db = Otter(otter_path)
        else:
            raise ValueError("Either the db or otter_path arguments must be provided!")

    res = db.query(**kwargs)

    if len(res) > result_length_tol:
        raise RuntimeError(
            f"This query returned {len(res)} results which is greater than the given "
            + f"tolerance of {result_length_tol}! Either increase the result_length_tol"
            + " keyword or pass in a stricter query!"
        )

    figs = []
    for t in res:
        try:
            fig = quick_view(
                t, ptype, sed_dim, dt_over_t, plotting_kwargs, **phot_cleaning_kwargs
            )
        except (KeyError, FailedQueryError):
            warn(f"No photometry associated with {t.default_name}, skipping!")
            continue

        fig.suptitle(t.default_name)
        figs.append(fig)

    return figs


def quick_view(
    t: Transient,
    ptype: str = "both",
    sed_dim: str = "freq",
    dt_over_t: float = 0,
    plotting_kwargs: dict = {},
    **kwargs,
) -> plt.Figure:
    """
    Generate a quick view (not necessarily publication ready) of the transients light
    curve, SED, or both. Default is to do both.

    Args:
        t (otter.Transient) : An otter Transient object to grab photometry from
        ptype (str) : The plot type to generate. Valid options are
                      - both -> Plot both light curve and sed (default)
                      - sed -> Plot just the sed
                      - lc -> Plot just the light curve
        sed_dim (str) : The x dimension to plot in the SED. Options are "freq" or
                        "wave". Default is "freq".
        dt_over_t (float) : The tolerance to split the days by. Default is 1 day. must
                           be unitless.
        plotting_kwargs (dict) : dictionary of key word arguments to pass to
                                 otter.plotter.plot_light_curve or
                                 otter.plotter.plot_sed.
        **kwargs : Any other arguments to pass to otter.Transient.clean_photometry

    Returns:
        The matplotlib figure used for plotting.
    """
    backend = plotting_kwargs.get("backend", "matplotlib.pyplot")
    if backend not in {"matplotlib.pyplot", "matplotlib"}:
        raise ValueError(
            "Only matplotlib.pyplot backend is available for quick_view!"
            + " To use plotly, use the plotting functionality individually!"
        )

    allphot = t.clean_photometry(**kwargs)
    allphot = allphot.sort_values("converted_date")
    allphot["time_tol"] = dt_over_t * allphot["converted_date"]
    allphot["time_diff"] = allphot["converted_date"].diff().fillna(-np.inf)
    allphot["time_grp"] = (allphot.time_diff > allphot.time_tol).cumsum()

    plt_lc = (ptype == "both") or (ptype == "lc")
    plt_sed = (ptype == "both") or (ptype == "sed")

    if ptype == "both":
        fig, (lc_ax, sed_ax) = plt.subplots(1, 2)
    elif ptype == "sed":
        fig, sed_ax = plt.subplots()
    elif ptype == "lc":
        fig, lc_ax = plt.subplots()

    if np.all(pd.isna(allphot.converted_flux_err)):
        flux_err = None
    else:
        flux_err = allphot.converted_flux_err

    if plt_lc:
        for filt, phot in allphot.groupby("filter_name"):
            plot_light_curve(
                date=phot.converted_date,
                flux=phot.converted_flux,
                flux_err=flux_err[allphot.filter_name == filt],
                xlabel=f"Date [{phot.converted_date_unit.values[0]}]",
                ylabel=f"Flux [{phot.converted_flux_unit.values[0]}]",
                ax=lc_ax,
                label=filt,
                **plotting_kwargs,
            )

    if plt_sed:
        for grp_name, phot in allphot.groupby("time_grp"):
            if sed_dim == "wave":
                wave_or_freq = phot.converted_wave
                xlab = f"Wavelength [{phot.converted_wave_unit.values[0]}]"
            elif sed_dim == "freq":
                wave_or_freq = phot.converted_freq
                xlab = f"Frequency [{phot.converted_freq_unit.values[0]}]"
            else:
                raise ValueError("sed_dim value is not recognized!")

            plot_sed(
                wave_or_freq=wave_or_freq,
                flux=phot.converted_flux,
                flux_err=flux_err[allphot.time_grp == grp_name],
                ax=sed_ax,
                xlabel=xlab,
                ylabel=f"Flux [{phot.converted_flux_unit.values[0]}]",
                label=phot.converted_date.mean(),
                **plotting_kwargs,
            )

        sed_ax.set_xscale("log")

    return fig


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
        fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel)

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
        fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel)

    return fig
