"""
Test the OTTER DataFinder class
"""

import pytest
import pandas as pd
from otter import DataFinder
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.utils import TableList
from requests.exceptions import ConnectionError, Timeout


def construct_data_finder(**kwargs):
    """
    Construct a DataFinder object and return it

    Based on the ASASSN-14li coordinates
    """
    ra = 192.06347
    dec = 17.77402
    return DataFinder(ra, dec, "deg", "deg", **kwargs)


def test_data_finder_constructor():
    """
    Test the constructor and make sure the instance variables are set correctly
    """
    # first, very basic, test
    df1 = construct_data_finder()
    assert isinstance(df1.coord, SkyCoord), "Not constructing a SkyCoord"
    assert df1.coord.ra.value == 192.06347, "Incorrect RA"
    assert df1.coord.dec.value == 17.77402, "Incorrect Dec"
    assert df1.name is None
    assert df1.z is None
    assert df1.redshift is None
    assert df1.bibcodes is None

    # next, more complex test
    df2 = construct_data_finder(
        name="ASASSN-14li",
        reference=["my_fake_reference"],
        redshift=0.1,  # this is made up for testing
    )
    assert isinstance(df1.coord, SkyCoord), "Not constructing a SkyCoord"
    assert df2.coord.ra.value == 192.06347, "Incorrect RA"
    assert df2.coord.dec.value == 17.77402, "Incorrect Dec"
    assert df2.name == "ASASSN-14li"
    assert df2.z == 0.1
    assert df2.redshift == 0.1
    assert df2.bibcodes == ["my_fake_reference"]


def test_repr():
    """
    Make sure the string representation is correct
    """
    df1 = construct_data_finder()
    assert str(df1) == "No Name DataFinder @ (RA, Dec)=(192.06347 deg,17.77402 deg)"

    df2 = construct_data_finder(
        name="ASASSN-14li",
        reference=["my_fake_reference"],
        redshift=0.1,  # this is made up for testing
    )
    assert str(df2) == "ASASSN-14li @ (RA, Dec)=(192.06347 deg,17.77402 deg)"


def test_iter():
    """
    Test the __iter__ method
    """
    keys1 = {"host_ra", "host_dec", "host_ra_units", "host_dec_units"}
    df1 = construct_data_finder()
    assert len(list({k for k, _ in df1} - keys1)) == 0

    keys2 = {
        "host_ra",
        "host_dec",
        "host_ra_units",
        "host_dec_units",
        "host_name",
        "host_redshift",
        "reference",
    }
    df2 = construct_data_finder(
        name="ASASSN-14li",
        reference=["my_fake_reference"],
        redshift=0.1,  # this is made up for testing
    )
    assert len(list({k for k, _ in df2} - keys2)) == 0


def test_query_simbad():
    """
    Test the query_simbad method
    """
    df1 = construct_data_finder()
    res = df1.query_simbad()

    assert len(res) == 3


def test_query_vizier():
    """
    Test the query_vizier method
    """
    df1 = construct_data_finder()
    res = df1.query_vizier()

    assert isinstance(res, TableList)
    assert len(res) >= 20, "missing some data"
    assert "II/294/sdss7" in res.keys()


def broken_test_query_atlas():
    """
    Test the query_atlas method

    THIS IS BROKEN ON GITHUB ACTIONS BECAUSE WE CANT SET THE ENV VARS!!!

    This one will inherently take a while!!! But it's because of atlas not us...
    """
    df1 = construct_data_finder()
    res = df1.query_atlas()

    assert "mjd" in res.columns
    assert "uJy" in res.columns
    assert "duJy" in res.columns
    assert len(res) >= 96, "missing some data!"


def test_query_ptf():
    """
    Test the method that queries palomer transient facility
    """
    df1 = construct_data_finder()
    res = df1.query_ptf()

    assert isinstance(res, Table)
    assert "obsmjd" in res.columns
    assert "mag_autocorr" in res.columns
    assert "meanmag" in res.columns
    assert len(res) >= 32, "Missing some PTF data!"


def test_query_ztf():
    """
    Test the method that queries the zwicky transient facility
    """
    df1 = construct_data_finder()
    res = df1.query_ztf()

    assert isinstance(res, Table)
    assert "mjd" in res.columns
    assert "mag" in res.columns
    assert "ra" in res.columns
    assert len(res) >= 1296, "Missing some ZTF data!"


def test_query_asassn():
    """
    Test the method that queries the ASASSN survey
    """
    df1 = construct_data_finder()
    res = df1.query_asassn()

    assert isinstance(res, pd.DataFrame)
    assert "asas_sn_id" in res.columns
    assert "jd" in res.columns
    assert "flux" in res.columns
    assert len(res) >= 1158, "Missing some ZTF data!"


def test_query_wise():
    """
    Test querying wise
    """
    df1 = construct_data_finder()
    res = df1.query_wise()

    assert isinstance(res, pd.DataFrame)
    assert "flux" in res.columns
    assert "date_mjd" in res.columns
    assert "name" in res.columns
    assert len(res) >= 35, "Missing some WISE data!"


@pytest.mark.skip(
    reason="There is a known issue with setuptools that makes this break\
when run in most conda environments!"
)
def test_query_alma():
    """
    Test querying ALMA
    """
    df1 = construct_data_finder()
    res = df1.query_alma()

    assert isinstance(res, Table)
    assert "obs_id" in res.columns
    assert "access_url" in res.columns


def test_query_first():
    """
    Test queryiung the FIRST radio survey
    """
    df1 = construct_data_finder()
    res = df1.query_first()

    assert isinstance(res, TableList)
    assert len(res) >= 1
    assert "J/ApJ/737/45/table2" in res.keys()


def test_query_nvss():
    """
    Test querying the NVSS radio survey
    """
    df1 = construct_data_finder()
    res = df1.query_nvss()

    assert isinstance(res, TableList)
    assert len(res) >= 1
    assert "J/ApJ/737/45/table1" in res.keys()


@pytest.mark.skip(
    reason="""This has a tendancy to timeout, not sure why but it sounds like its
    not on us (hopefully)..."""
)
def test_query_sparcl():
    """
    Test querying SPARCL for spectra

    Note: Sparcl likes to throttle us so I'm wrapping this in a try/except for some very
    specific errors that requests will throw
    """

    df1 = construct_data_finder()
    try:
        res = df1.query_sparcl()
        assert isinstance(res, Table)
        assert len(res) >= 1, "Missing some spectroscopic data!"

    except (ConnectionError, Timeout):  # NOTE: NEVER add AttributeError to this list
        pass


def test_query_heasarc():
    """
    Test querying the Heasarc x-ray catalog
    """

    df1 = construct_data_finder()

    # test with x-ray
    res = df1.query_heasarc(catalog="xray")
    assert isinstance(res, Table)
    assert len(res) >= 14, "Missing some HEASARC data"

    # test with radio
    res2 = df1.query_heasarc(catalog="radio")
    assert isinstance(res2, Table)
    assert len(res2) >= 2, "Missing some HEASARC data"
