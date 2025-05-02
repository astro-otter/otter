"""
Unit tests for the otter.Otter class methods
"""

import os
from otter import Otter, Transient
from otter.exceptions import FailedQueryError
from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
import pandas as pd
import pytest

OTTER_URL = os.environ.get("OTTER_TEST_URL")
OTTER_TEST_PASSWORD = os.environ.get("OTTER_TEST_PASSWORD")


def test_otter_constructor():
    """
    Just make sure everything constructs correctly
    """

    db = Otter(url=OTTER_URL, password=OTTER_TEST_PASSWORD)
    assert isinstance(db, Otter)


def test_get_meta():
    """
    Tests the Otter.get_meta method and make sure it returns as expected
    """

    db = Otter(url=OTTER_URL, password=OTTER_TEST_PASSWORD)

    # first make sure everything is just copied over correctly
    allmeta = db.get_meta()
    true_keys = ["name", "coordinate", "date_reference", "distance", "classification"]

    assert all(k in d for d in allmeta for k in true_keys)

    # now we can try real queries
    metahyz = db.get_meta(names="2018hyz")[0]
    assert isinstance(metahyz, Transient)
    assert metahyz["name/default_name"] == "2018hyz"
    assert metahyz["date_reference"][0]["value"] == "2018-10-14"
    assert metahyz["date_reference"][0]["date_format"] == "iso"
    assert metahyz["classification"][0]["object_class"] == "TDE"


def test_cone_search():
    """
    Tests the Otter.cone_search method
    """

    db = Otter(url=OTTER_URL, password=OTTER_TEST_PASSWORD)

    # just search around '2018hyz' coordinates to make sure it picks it up
    coord = SkyCoord(151.711964138, 1.69279894089, unit="deg")
    res = db.cone_search(coord)[0]
    assert res["name/default_name"] == "2018hyz"


def test_get_phot():
    """
    Tests the Otter.get_phot method

    We know from the transients.clean_photometry tests that the conversions
    work as expected. So, this will just test that everything comes out as expected.
    """

    db = Otter(url=OTTER_URL, password=OTTER_TEST_PASSWORD)

    true_keys = [
        "name",
        "converted_flux",
        "converted_flux_err",
        "converted_date",
        "converted_wave",
        "converted_freq",
        "converted_flux_unit",
        "converted_date_unit",
        "converted_wave_unit",
        "converted_freq_unit",
        "obs_type",
        "upperlimit",
        "reference",
    ]

    names = ["2018hyz", "2018zr", "ASASSN-14li"]

    # first with returning an astropy table (the default)
    allphot = db.get_phot(names=names)
    assert isinstance(allphot, Table)
    assert all(k in allphot.keys() for k in true_keys)
    assert len(np.unique(allphot["converted_flux_unit"])) == 1
    assert allphot["converted_flux_unit"][0] == "mag(AB)"

    # then with returning a pandas DataFrame
    allphot = db.get_phot(names=names, return_type="pandas")
    assert isinstance(allphot, pd.DataFrame)
    assert all(k in allphot for k in true_keys)
    assert len(np.unique(allphot.converted_flux_unit)) == 1
    assert allphot.converted_flux_unit.iloc[0] == "mag(AB)"

    # then make sure it throws the FailedQueryError
    with pytest.raises(FailedQueryError):
        db.get_phot(names="foo")


def test_query():
    """
    Tests the Otter.query method that basically all of this is based on

    A lot of these have been tested in other unit tests in thie file
    but lets make sure it's complete
    """

    db = Otter(url=OTTER_URL, password=OTTER_TEST_PASSWORD)

    # test min and max z queries
    zgtr1 = db.query(minz=1)
    assert len(zgtr1) >= 2
    true_result = ["Sw J2058+05", "2022cmc", "CXOU J0332"]
    assert all(t["name/default_name"] in true_result for t in zgtr1)

    zless001 = db.query(maxz=0.001)
    result = ["NGC 247", "IGR J17361-4441"]
    assert all(t["name/default_name"] in result for t in zless001)

    # test refs
    # res = db.query(refs="2020MNRAS.tmp.2047S")[0]
    # assert res["name/default_name"] == "2018hyz"

    # test hasphot and hasspec
    assert len(db.query(hasspec=True)) == 0
    assert "ASASSN-20il" not in {t["name/default_name"] for t in db.query(hasphot=True)}
