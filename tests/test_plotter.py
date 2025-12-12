"""
Some testing for the mini plotting module that's part of this api

Author: Noah Franz
"""

import matplotlib.pyplot as plt
import plotly
from otter import Otter
from otter.plotter import plotter
import pytest


def _get_test_data():
    db = Otter()
    t = db.query(names="ASASSN-14li")
    return t[0]


def test_otter_plotter_class():
    # try initializing the class a bunch of ways
    op = plotter.OtterPlotter(backend="matplotlib.pyplot")
    assert op.backend == "matplotlib.pyplot"

    op = plotter.OtterPlotter(backend="plotly.graph_objects")
    assert op.backend == "plotly.graph_objects"

    op = plotter.OtterPlotter(backend="plotly")
    assert op.backend == "plotly.graph_objects"

    op = plotter.OtterPlotter(backend="matplotlib")
    assert op.backend == "matplotlib.pyplot"

    with pytest.raises(ValueError):
        plotter.OtterPlotter(backend="test")

    # then test the two backends with some basics
    op = plotter.OtterPlotter(backend="matplotlib")
    ax = op.plot([1], [1])
    assert isinstance(ax, plt.Axes)

    op = plotter.OtterPlotter(backend="plotly")
    fig = op.plot([1], [1])
    assert isinstance(fig, plotly.graph_objs._figure.Figure)


def test_plot_light_curve():
    t = _get_test_data()
    phot = t.clean_photometry(obs_type="uvoir", flux_unit="mag(AB)")

    # first try with matplotlib
    ax = plotter.plot_light_curve(
        phot.converted_date,
        phot.converted_flux,
        date_err=None,
        flux_err=phot.converted_flux_err,
        backend="matplotlib",
    )
    assert isinstance(ax, plt.Axes)

    # then try with plotly
    fig = plotter.plot_light_curve(
        phot.converted_date,
        phot.converted_flux,
        date_err=None,
        flux_err=phot.converted_flux_err,
        backend="plotly",
    )
    assert isinstance(fig, plotly.graph_objs._figure.Figure)


def test_plot_sed():
    t = _get_test_data()
    phot = t.clean_photometry(obs_type="uvoir", flux_unit="mag(AB)")

    ax = plotter.plot_sed(
        phot.converted_wave,
        phot.converted_flux,
        flux_err=phot.converted_flux_err,
        backend="matplotlib",
        ylabel="Flux [AB]",
        xlabel="Wavelength [nm]",
    )
    assert isinstance(ax, plt.Axes)
    assert ax.get_ylabel() == "Flux [AB]"
    assert ax.get_xlabel() == "Wavelength [nm]"

    fig = plotter.plot_sed(
        phot.converted_freq,
        phot.converted_flux,
        flux_err=phot.converted_flux_err,
        backend="plotly",
        ylabel="Flux [AB]",
        xlabel="Frequency [GHz]",
    )
    assert isinstance(fig, plotly.graph_objs._figure.Figure)
    assert fig.layout.xaxis.title.text == "Frequency [GHz]"
    assert fig.layout.yaxis.title.text == "Flux [AB]"


def test_quick_view():
    t = _get_test_data()

    fig = plotter.quick_view(t, ptype="both", sed_dim="wave", obs_type="uvoir")

    assert isinstance(fig, plt.Figure)
    assert len(fig.get_axes()) == 2
    lc_ax, sed_ax = fig.get_axes()
    assert lc_ax.get_xlabel() == "Date [MJD]"
    assert lc_ax.get_ylabel() == "Flux [mag(AB)]"
    assert sed_ax.get_xlabel() == "Wavelength [nm]"
    assert sed_ax.get_ylabel() == "Flux [mag(AB)]"

    fig = plotter.quick_view(t, ptype="lc", obs_type="radio", flux_unit="mJy")
    assert isinstance(fig, plt.Figure)
    assert len(fig.get_axes()) == 1
    ax = fig.get_axes()[0]
    assert ax.get_xlabel() == "Date [MJD]"
    assert ax.get_ylabel() == "Flux [mJy]"

    fig = plotter.quick_view(t, ptype="sed", obs_type="uvoir")
    assert isinstance(fig, plt.Figure)
    assert len(fig.get_axes()) == 1
    ax = fig.get_axes()[0]
    assert ax.get_xlabel() == "Frequency [GHz]"
    assert ax.get_ylabel() == "Flux [mag(AB)]"

    with pytest.raises(ValueError):
        plotter.quick_view(t, plotting_kwargs=dict(backend="plotly"))


def test_query_quick_view():
    db = Otter()

    fig = plotter.query_quick_view(
        db,
        ptype="both",
        sed_dim="wave",
        phot_cleaning_kwargs=dict(obs_type="uvoir"),
        names="ASASSN-14li",
    )

    assert all(isinstance(f, plt.Figure) for f in fig)
    assert len(fig[0].get_axes()) == 2
    lc_ax, sed_ax = fig[0].get_axes()
    assert lc_ax.get_xlabel() == "Date [MJD]"
    assert lc_ax.get_ylabel() == "Flux [mag(AB)]"
    assert sed_ax.get_xlabel() == "Wavelength [nm]"
    assert sed_ax.get_ylabel() == "Flux [mag(AB)]"

    fig = plotter.query_quick_view(
        db,
        ptype="lc",
        phot_cleaning_kwargs=dict(obs_type="radio", flux_unit="mJy"),
        names="Sw J1644+57",
    )
    assert isinstance(fig[0], plt.Figure)
    assert len(fig[0].get_axes()) == 1
    ax = fig[0].get_axes()[0]
    assert ax.get_xlabel() == "Date [MJD]"
    assert ax.get_ylabel() == "Flux [mJy]"

    fig = plotter.query_quick_view(
        db,
        ptype="sed",
        phot_cleaning_kwargs=dict(obs_type="uvoir"),
        names="ASASSN-14li",
    )
    assert isinstance(fig[0], plt.Figure)
    assert len(fig[0].get_axes()) == 1
    ax = fig[0].get_axes()[0]
    assert ax.get_xlabel() == "Frequency [GHz]"
    assert ax.get_ylabel() == "Flux [mag(AB)]"

    with pytest.raises(RuntimeError):
        fig = plotter.query_quick_view(db, minz=1, maxz=2, result_length_tol=1)
