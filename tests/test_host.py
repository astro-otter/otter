"""
Test the otter Host class
"""

from otter import Host
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.utils import TableList


def test_constructor():
    """
    Test the Host constructor
    """
    host_ra = 192.06249999999997
    host_dec = 17.77401111111111
    host_name = "SDSS J124815.23+174626.4"
    h = Host(
        host_ra=host_ra,
        host_dec=host_dec,
        host_name=host_name,
        host_ra_units="deg",
        host_dec_units="deg",
        transient_name="ASASSN-14li",
        reference=["foobar"],
        host_redshift=0.1,
    )

    assert isinstance(h.coord, SkyCoord), "Not constructing a SkyCoord"
    assert h.coord.ra.value == host_ra, "Incorrect RA"
    assert h.coord.dec.value == host_dec, "Incorrect Dec"
    assert h.name == host_name
    assert h.z == 0.1
    assert h.redshift == 0.1
    assert h.bibcodes == ["foobar"]
    assert h.transient_name == "ASASSN-14li"


def test_pcc(fake_astroquery):
    """
    Test the PCC score testing (the VizieR SDSS magnitude lookup is mocked)
    """
    fake_astroquery.result = TableList(
        {
            "II/294/sdss7": Table({"rmag": [20.0]}),
            "V/147/sdss12": Table({"rmag": [16.6]}),  # should use the latest SDSS
        }
    )

    host_ra = 192.06249999999997
    host_dec = 17.77401111111111
    host_name = "SDSS J124815.23+174626.4"
    h = Host(
        host_ra=host_ra,
        host_dec=host_dec,
        host_name=host_name,
        host_ra_units="deg",
        host_dec_units="deg",
        transient_name="ASASSN-14li",
        reference=["foobar"],
        host_redshift=0.1,
    )

    # get the testing path
    coord = SkyCoord(192.06343875, 17.77402083, unit="deg")
    pcc = h.pcc(coord)

    assert pcc < 0.001, "pcc seems suspiciously high for ASASSN-14li as a test case..."

    pcc2 = h.pcc(SkyCoord(0, 0, unit="deg"))
    assert pcc2 > 0.9, "pcc seems too small for this input coordinate"


def test_query_blast(fake_blast):
    """
    Test the query of BLAST (the BLAST API is mocked)
    """

    host = Host.query_blast(tns_name="AT2024advm")

    assert "name=2024advm" in fake_blast[0]
    assert host.transient_name == "2024advm", "name incorrect!"
    assert host.name == "FAKE-HOST-1"
    assert host.z == 0.05, "should prefer the spectroscopic redshift"
