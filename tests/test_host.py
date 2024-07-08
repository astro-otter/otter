"""
Test the otter Host class
"""

import os
from otter import Host
from otter import Otter
from astropy.coordinates import SkyCoord


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


def test_pcc():
    """
    Test the PCC score testing
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

    # get the testing path
    otterpath = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), ".otter-testing"
    )
    db = Otter(otterpath)
    t = db.query(names=h.transient_name)[0]
    pcc = h.pcc(t.get_skycoord())

    assert pcc < 0.001, "pcc seems suspiciously high for ASASSN-14li as a test case..."

    pcc2 = h.pcc(SkyCoord(0, 0, unit="deg"))
    assert pcc2 > 0.9, "pcc seems too small for this input coordinate"
