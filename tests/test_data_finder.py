"""
Test the OTTER DataFinder class
"""

import os
import pytest
import shutil
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


def test_query_simbad(fake_astroquery):
    """
    Test the query_simbad method (SIMBAD is mocked)
    """
    from astroquery.simbad import Simbad

    fake_astroquery.result = Table({"main_id": ["ASASSN-14li"]})
    df1 = construct_data_finder()
    res = df1.query_simbad()

    assert len(res) == 1
    module, args, kwargs = fake_astroquery.calls[0]
    assert module is Simbad
    assert args[0] == df1.coord
    assert kwargs["radius"] == "5 arcsec"


def test_query_vizier(fake_astroquery):
    """
    Test the query_vizier method (VizieR is mocked)
    """
    from astroquery.vizier import Vizier
    from otter.io.data_finder import VIZIER_LARGE_CATALOGS

    fake_astroquery.result = TableList({"II/294/sdss7": Table({"rmag": [17.0]})})
    df1 = construct_data_finder()
    res = df1.query_vizier()

    assert isinstance(res, TableList)
    assert "II/294/sdss7" in res.keys()
    module, _, kwargs = fake_astroquery.calls[0]
    assert module is Vizier
    assert kwargs["catalog"] == VIZIER_LARGE_CATALOGS


class _FakeATLASResponse:
    def __init__(self, status_code, json_data=None, text=""):
        self.status_code = status_code
        self._json = json_data
        self.text = text

    def json(self):
        return self._json


class _FakeATLASSession:
    """Scripted ATLAS forced photometry server: queue -> finished task -> result"""

    task_url = "https://fallingstar-data.com/forcedphot/queue/1/"
    result_url = "https://fallingstar-data.com/forcedphot/static/results/1.txt"
    result_text = "\n".join(
        [
            "###MJD m dm uJy duJy F err chi/N RA Dec x y maj min phi apfit "
            "mag5sig Sky Obs"
        ]
        + [
            f"{60000 + i:.5f} 18.0 0.05 230.0 15.0 {f} 0 1.2 192.06 17.77 "
            "5000 5000 2.0 2.0 0 -0.4 19.5 20.0 01a"
            for i in range(4)
            for f in "co"
        ]
    )

    def __init__(self):
        self.posts = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def post(self, url, headers=None, data=None):
        self.posts.append((url, headers, data))
        return _FakeATLASResponse(201, {"url": self.task_url})

    def get(self, url, headers=None):
        if url == self.task_url:
            return _FakeATLASResponse(
                200,
                {"finishtimestamp": "now", "result_url": self.result_url},
            )
        assert url == self.result_url
        return _FakeATLASResponse(200, text=self.result_text)


def test_query_atlas(monkeypatch):
    """
    Test the query_atlas method (the ATLAS forced photometry server is mocked)
    """
    session = _FakeATLASSession()
    monkeypatch.setenv("ATLAS_API_TOKEN", "fake-token")
    monkeypatch.setattr("otter.io.data_finder.requests.Session", lambda: session)

    df1 = construct_data_finder()
    res = df1.query_atlas(days_ago=7)

    assert "mjd" in res.columns
    assert "uJy" in res.columns
    assert "duJy" in res.columns
    assert len(res) > 0
    assert session.posts[0][1]["Authorization"] == "Token fake-token"
    assert session.posts[0][2]["ra"] == 192.06347
    assert session.posts[0][2]["dec"] == 17.77402


def test_query_ptf(fake_astroquery):
    """
    Test the method that queries palomer transient facility (IRSA is mocked)
    """
    from astroquery.ipac.irsa import Irsa

    fake_astroquery.result = Table(
        {"obsmjd": [55000.0], "mag_autocorr": [18.0], "meanmag": [18.1]}
    )
    df1 = construct_data_finder()
    res = df1.query_ptf()

    assert isinstance(res, Table)
    assert "obsmjd" in res.columns
    assert "mag_autocorr" in res.columns
    assert "meanmag" in res.columns
    module, _, kwargs = fake_astroquery.calls[0]
    assert module is Irsa
    assert kwargs["catalog"] == "ptf_lightcurves"


def test_query_ztf(monkeypatch):
    """
    Test the method that queries the zwicky transient facility (IRSA is mocked)
    """
    import io
    from astropy.io.votable import from_table
    from conftest import FakeHTTPResponse

    buf = io.BytesIO()
    from_table(
        Table({"mjd": [58000.0, 58001.0], "mag": [18.0, 18.1], "ra": [192.06, 192.06]})
    ).to_xml(buf)

    urls = []

    def fake_urlopen(url):
        urls.append(url)
        return FakeHTTPResponse(buf.getvalue())

    monkeypatch.setattr("otter.io.data_finder.urlopen", fake_urlopen)

    df1 = construct_data_finder()
    res = df1.query_ztf()

    assert isinstance(res, Table)
    assert "mjd" in res.columns
    assert "mag" in res.columns
    assert "ra" in res.columns
    assert len(res) == 2
    assert "POS=CIRCLE%20192.06347%2017.77402" in urls[0]


def test_query_asassn(monkeypatch):
    """
    Test the method that queries the ASASSN survey (SkyPatrol is mocked)
    """
    calls = []

    class FakeSkyPatrolClient:
        def cone_search(self, ra, dec, **kwargs):
            calls.append((ra, dec, kwargs))

            class LightCurve:
                data = pd.DataFrame(
                    {
                        "asas_sn_id": [1, 1],
                        "jd": [2458000.5, 2458001.5],
                        "flux": [1.0, 2.0],
                    }
                )

            return LightCurve()

    monkeypatch.setattr("pyasassn.client.SkyPatrolClient", FakeSkyPatrolClient)

    df1 = construct_data_finder()
    res = df1.query_asassn()
    assert isinstance(res, pd.DataFrame)
    assert "asas_sn_id" in res.columns
    assert "jd" in res.columns
    assert "flux" in res.columns
    assert calls[0][:2] == (192.06347, 17.77402)
    assert calls[0][2]["radius"] == 5.0


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


def test_query_first(fake_astroquery):
    """
    Test queryiung the FIRST radio survey (VizieR is mocked)
    """
    fake_astroquery.result = TableList({"J/ApJ/737/45/table2": Table({"Fp": [1.0]})})
    df1 = construct_data_finder()
    res = df1.query_first()

    assert isinstance(res, TableList)
    assert "J/ApJ/737/45/table2" in res.keys()
    assert fake_astroquery.calls[0][2]["catalog"] == "J/ApJ/737/45/table2"

    # res, res_img = df1.query_first(get_image=True)
    # assert isinstance(res_img,


def test_query_nvss(fake_astroquery):
    """
    Test querying the NVSS radio survey (VizieR is mocked)
    """
    fake_astroquery.result = TableList({"J/ApJ/737/45/table1": Table({"S1.4": [1.0]})})
    df1 = construct_data_finder()
    res = df1.query_nvss()

    assert isinstance(res, TableList)
    assert "J/ApJ/737/45/table1" in res.keys()
    assert fake_astroquery.calls[0][2]["catalog"] == "J/ApJ/737/45/table1"


@pytest.mark.skip(
    reason="DataLab & Sparcl have a tendency to timeout, which isn't on us!"
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


@pytest.mark.skip(reason="pyvo is rejecting our queries of the heasarc catalog")
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


@pytest.mark.skip(
    "This works, but fails on MacOS *only* on github actions, for some reason..."
)
def test_query_wiserep():
    """
    Test querying WISeREP for spectra
    """

    df = construct_data_finder()
    with pytest.raises(ValueError):
        df.query_wiserep(iau_name="ASASSN-14li")

    df1 = DataFinder(ra=0, dec=0, ra_units="deg", dec_units="deg", name="2022cmc")

    # then test with 22cmc since we know that works well
    meta, res = df1.query_wiserep()
    assert isinstance(meta, pd.DataFrame)
    assert isinstance(res, dict)
    assert isinstance(res[0], pd.DataFrame)
    assert not os.path.exists("spectra")

    meta, res = df1.query_wiserep(iau_name="AT2022cmc", rm_tmp_dir=False)
    assert isinstance(meta, pd.DataFrame)
    assert isinstance(res, dict)
    assert isinstance(res[0], pd.DataFrame)
    assert os.path.exists("spectra")

    # then cleanup
    shutil.rmtree("spectra")
