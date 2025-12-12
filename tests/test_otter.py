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


def test_otter_constructor():
    """
    Just make sure everything constructs correctly
    """

    db = Otter()
    assert isinstance(db, Otter)

    # and using a custom local directory
    otter_test_dir = os.path.join(os.getcwd(), "otter-test")
    db = Otter(
        datadir=otter_test_dir,
        gen_summary=True,
    )

    assert os.path.exists(otter_test_dir)  # make sure it was created
    assert isinstance(db, Otter)


def test_get_meta():
    """
    Tests the Otter.get_meta method and make sure it returns as expected
    """

    db = Otter()

    # first make sure everything is just copied over correctly
    allmeta = db.get_meta()
    req_keys = ["name", "coordinate"]

    assert all(k in d for d in allmeta for k in req_keys)

    # now we can try real queries
    true_keys = ["name", "coordinate", "date_reference", "distance", "classification"]
    metahyz = db.get_meta(names="2018hyz")[0]
    assert isinstance(metahyz, Transient)
    assert all(k in metahyz for k in true_keys)
    assert metahyz["name/default_name"] == "2018hyz"
    assert metahyz["date_reference"][0]["value"] == "2018-10-14"
    assert metahyz["date_reference"][0]["date_format"] == "iso"
    assert metahyz["classification"]["value"][0]["object_class"] == "TDE"


def test_cone_search():
    """
    Tests the Otter.cone_search method
    """

    db = Otter()

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

    db = Otter()

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

    names = ["2018hyz", "2018zr", "ASASSN-14li", "2016fnl", "J123715"]

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

    # some other random tests
    with pytest.raises(OSError):
        db.get_phot(names="ASASSN-14li", return_type="foobar")


def test_query():
    """
    Tests the Otter.query method that basically all of this is based on

    A lot of these have been tested in other unit tests in thie file
    but lets make sure it's complete
    """

    db = Otter()

    # test min and max z queries
    zgtr1 = db.query(minz=1)
    assert len(zgtr1) >= 2
    true_result = ["Swift J2058.4+0516", "2022cmc", "CXOU J0332"]
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

    # test has_*_phot
    assert len(db.query(has_radio_phot=True)) >= 92
    assert len(db.query(has_xray_phot=True)) >= 35
    assert len(db.query(has_uvoir_phot=True)) >= 120
    assert len(db.query(has_radio_phot=True, has_xray_phot=True)) < len(
        db.query(has_radio_phot=True)
    )

    # test classification related queries
    assert len(db.query(spec_classed=True)) > 140
    assert len(db.query(unambiguous=True)) > 190
    assert len(db.query(classification="SLSN")) >= 1

    # check that querying based on references works
    fake_test_bibcodes = ["'tasdfasdf...lkjsfd'", "'asfd...kjasdf...lkjs'"]
    assert len(db.query(refs=fake_test_bibcodes)) == 0
    assert len(db.query(refs=fake_test_bibcodes[0])) == 0

    # try querying with the fake "private" data
    db2 = Otter(
        datadir=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "private_otter_data"
        ),
        gen_summary=True,
    )
    res = db2.query(names="2018hyz", query_private=True)
    assert res[0].default_name == "2018hyz"

    # check that some errors are thrown when appropriate
    with pytest.raises(Exception):
        db.query(names=10)
        db.query(refs=10)

    # try just retrieving all of the private data
    all_res = db2._query_datadir()
    assert len(all_res) == 1

    # check some other basic queries agains the private data
    res = db2.query(names=["2018hyz", "2022cmc"])
    assert len(res) == 2

    with pytest.raises(AttributeError):
        res = db2.query(names="2018hyz", hasspec=True, query_private=True)

    res = db2.query(names="2018hyz", hasphot=True, query_private=True)
    assert len(res) == 1

    res = db2.query(names="2018hyz", maxz=1, query_private=True)
    assert len(res) == 1

    res = db2.query(minz=1, query_private=True)
    assert len(res) >= 3
    assert "2018hyz" not in [t.default_name for t in res]

    res = db2._query_datadir(refs="2018TNSCR1764....1A")
    assert len(res) == 1

    res = db2._query_datadir(refs=["2018TNSCR1764....1A", "2018TNSTR1708....1B"])
    assert len(res) == 1

    # THESE ALL MUST HAPPEN AT THE END OF THE QUERYING
    # they are testing what happens if we mess with the summary table in the
    # private data directory
    summary_table_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "private_otter_data", "summary.csv"
    )
    summary_table = pd.read_csv(summary_table_path)
    summary_table.head(0).to_csv(summary_table_path)
    res = db2._query_datadir(names="2018hyz")
    assert isinstance(res, list)
    assert len(res) == 0

    # now delete the summary file because it should auto regenerate
    os.remove(summary_table_path)
    res = db2._query_datadir(names="2018hyz")
    assert len(res) == 1

    # test passing in lists and arrays of names
    res = db.query(names=["2018hyz", "2019azh"])
    assert len(res) == 2

    res = db.query(names=np.array(["2018hyz", "2019azh"]))
    assert len(res) == 2


def test_from_csvs():
    """
    This tests the "from_csvs" method which allows for interaction with locally
    stored datasets
    """

    metapath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meta-test.csv")
    photpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "phot-test.csv")
    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "from-csvs-out")

    db = Otter.from_csvs(metafile=metapath, photfile=photpath, local_outpath=outpath)

    assert isinstance(db, Otter)

    metapath = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "gaia16aax-meta.csv"
    )
    photpath = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "gaia16aax-phot.csv"
    )
    outpath = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "from-csvs-out-gaia16aax"
    )

    db = Otter.from_csvs(metafile=metapath, photfile=photpath, local_outpath=outpath)

    assert isinstance(db, Otter)
