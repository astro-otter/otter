"""
Test the transient getter methods. I distinguish this from some of the other
methods (like the __add__ and corresponding private methods) since they
are more complex and require additional tests.
"""

import pytest
import numpy as np
from otter import Transient, Host
from otter.exceptions import OtterLimitationError, IOError
from astropy.coordinates import SkyCoord
from astropy.time import Time


def test_transient_constructor():
    """
    Test the constructor and make sure we have the expected instance variables
    """
    test_json = generate_test_json()

    # make sure the constructor works
    t = Transient(test_json)

    # test that the default_name attribute is correct
    assert t.default_name == "Sw J1644+57", "default_name is incorrect!"

    # test that the srcmap is correct
    # choose a few at random
    srcmap_msg = "Something with the srcmap is wrong!"
    assert t.srcmap["TNS"] == "TNS", srcmap_msg
    assert t.srcmap["2012MNRAS.421.1942W"] == "Wiersema et al. (2012)", srcmap_msg
    assert t.srcmap["2011ApJ...737..103S"] == "Schlafly & Finkbeiner (2011)", srcmap_msg


def test_getitem():
    """
    Test the __getitem__ method
    """
    msg = "Something is broken with __getitem__"

    # create a transient object to test
    t = Transient(generate_test_json())

    # test the "classic" dictionary syntax
    assert t["name"]["default_name"] == "Sw J1644+57", msg
    assert t["name"]["alias"][-1]["value"] == "Swift J1644+57", msg
    assert t["coordinate"][0]["ra"] == "16 44 49.93130", msg

    # use the same tests but with the "/" syntax
    assert t["name/default_name"] == "Sw J1644+57", msg
    assert t["name/alias"][-1]["value"] == "Swift J1644+57", msg
    assert t["coordinate"][0]["ra"] == "16 44 49.93130", msg

    # test it where we pass in a list or tuple
    test_keys = ["name", "photometry"]
    out = t[test_keys]
    assert out["name/default_name"] == "Sw J1644+57", msg
    assert out["photometry"][0]["raw"][0] == 0.25, msg


def test_setitem():
    """
    Test the __setitem__ method
    """
    msg = "somethings wrong with setitem!"

    t = Transient()

    # first test is basic
    test1 = {"mytest": "nothing to see here"}
    t["test1"] = test1

    assert t["test1"] == test1, msg
    assert t["test1"]["mytest"] == "nothing to see here", msg

    # next test tries to set an item with the / syntax
    t["test1/myothertest"] = "this is cool!"
    assert t["test1"]["myothertest"] == "this is cool!", msg
    assert t["test1/myothertest"] == "this is cool!", msg

    # now try overwriting mytest
    t["test1/mytest"] = "this is better"
    assert t["test1"]["mytest"] == "this is better", msg
    assert t["test1/mytest"] == "this is better", msg


def test_delitem():
    """
    Test the __delitem__ method
    """

    msg = "something went wrong with __delitem__ method!"

    # test with the classical syntax
    t = Transient({"test1": {"mytest": "nothing to see here"}})

    del t["test1"]["mytest"]
    assert t["test1"] == {}, msg

    # test with the '/' syntax
    # test with the classical syntax
    t = Transient({"test1": {"mytest": "nothing to see here"}})

    with pytest.raises(OtterLimitationError):
        del t["test1/mytest"]


def test_len():
    msg = "Something is wrong with the __len__ method!"

    t = Transient({"test1": {"mytest": "nothing to see here", "other": 1}})

    assert len(t) == 1, msg
    assert len(t["test1"]) == 2, msg


def test_repr():
    """
    Test the string representation of a Transient object
    """
    msg = "Something with the string representation is unexpected!"
    t = Transient(generate_test_json())

    str_repr = str(t)
    corr_str_repr = f"Transient(\n\tName: Sw J1644+57,\n\tKeys: {t.keys()}\n)"

    assert str_repr == corr_str_repr, msg


def test_iter():
    """
    Test that my __iter__ overwrite works
    """
    true_vals = [
        {"value": "Sw J1644+57", "reference": ["Swift"]},
        {"value": "GRB 110328A", "reference": [["2011Sci...333..203B"]]},
        {"value": "Swift J164449.3+573451", "reference": [["2011Sci...333..203B"]]},
        {
            "value": "Swift J1644+57",
            "reference": [["2017ApJ...838..149A", "2011Sci...333..203B"]],
        },
    ]
    to_iterate = Transient(generate_test_json())["name/alias"]
    for found, true in zip(to_iterate, true_vals):
        assert found == true, "iterating is not working!"


def test_keys():
    """
    Test the keys method
    """
    msg = "Something is wrong with the keys() method!"
    true_keys = [
        "schema_version",
        "distance",
        "filter_alias",
        "reference_alias",
        "date_reference",
        "name",
        "photometry",
        "coordinate",
        "classification",
        "host",
    ]

    t = Transient(generate_test_json())

    assert list(t.keys()) == true_keys, msg


def test_get_meta():
    """
    Test the get meta method of Transient
    """

    msg = "get_meta method is broken!"

    t = Transient(generate_test_json())
    meta = t.get_meta()
    meta2 = t.get_meta(keys=["name", "classification"])

    # check the keys are correct
    true_meta_keys = [
        "schema_version",
        "distance",
        "filter_alias",
        "reference_alias",
        "date_reference",
        "name",
        "coordinate",
        "classification",
        "host",
    ]
    assert list(meta.keys()) == true_meta_keys, msg
    assert list(meta2.keys()) == ["name", "classification"], msg

    # check that most of the values are correct
    # name first
    assert meta["name/default_name"] == "Sw J1644+57", msg
    assert meta["name/default_name"] == t.default_name, msg
    assert meta2["name/default_name"] == "Sw J1644+57", msg

    # now classification
    assert meta["classification"][0]["object_class"] == "TDE", msg
    assert meta2["classification"][0]["object_class"] == "TDE", msg

    # now distance
    assert meta["distance"][0]["value"] == "0.354", msg
    with pytest.raises(KeyError):
        meta2["distance"]  # make sure this key isn't in this metadata


def test_get_skycoord():
    """
    Test the get skycoord method
    """

    t = Transient(generate_test_json())

    skycoord = t.get_skycoord()

    assert isinstance(skycoord, SkyCoord), "get_skycoord did not return a SkyCoord!"
    assert str(skycoord.ra) == "251d12m28.9695s", "RA does not match!"
    assert str(skycoord.dec) == "57d34m59.6893s", "Dec does not match!"


def test_get_discovery_date():
    """
    test the get_discovery_date method
    """

    t = Transient(generate_test_json())
    discdate = t.get_discovery_date()

    assert isinstance(discdate, Time), "get_discovery_date did not return a Time!"
    assert str(discdate) == "2011-03-29 00:00:00.000", "The time string is incorrect!"


def test_get_redshift():
    """
    Test the get redshift method
    """

    t = Transient(generate_test_json())
    z = t.get_redshift()

    assert float(z) == 0.354, "get_redshift did not return the correct value!"


def test_get_host():
    """
    Test the get host method

    This is pretty bare bones because the rest of the output is tested with test_host
    """

    t = Transient(generate_test_json())

    host = t.get_host()
    assert isinstance(host, list)
    assert isinstance(host[0], Host)


def test_clean_photometry():
    """
    Make sure we can correctly clean the photometry
    """

    msg = "Something is broken with the photometry cleaning!"
    t = Transient(generate_test_json())

    # first with just the default options
    phot = t.clean_photometry()
    uq_obs_types = phot["obs_type"].unique()

    assert all(ot in uq_obs_types for ot in ["radio", "uvoir", "xray"]), msg
    assert len(phot["converted_flux_unit"].unique()) == 1, msg
    assert phot["converted_flux_unit"][0] == "mag(AB)", msg

    uvoir = phot[phot.obs_type == "uvoir"]
    assert np.isclose(phot["converted_flux"].iloc[0], 17.905, atol=1e-2), msg
    assert 17.59 in list(uvoir["converted_flux"]), msg

    # then with different "by"
    with pytest.raises(IOError):
        t.clean_photometry(by="foo")

    phot_val = t.clean_photometry(by="value").reset_index()

    assert len(phot_val.obs_type.unique()) == 1, msg
    assert phot_val.obs_type.iloc[0] == "xray", msg
    assert len(phot_val["converted_flux_unit"].unique()) == 1, msg
    assert phot_val["converted_flux_unit"][0] == "mag(AB)", msg

    # then with different returned units
    # and only get the radio data
    phot_non_default = t.clean_photometry(
        flux_unit="uJy", freq_unit="MHz", obs_type="radio"
    ).reset_index()

    assert len(phot_non_default.obs_type.unique()) == 1, msg
    assert phot_non_default.obs_type.iloc[0] == "radio", msg
    assert np.isclose(phot_non_default.converted_flux.iloc[0], 0.25e3), msg


# a test json file
# I use this here instead of reading in existing files
# because I want to be able to control exactly what goes
# into the test cases!
def generate_test_json():
    """
    Generate a test json file

    This is a mixture of data on Sw J1644 plus some random
    other photometry points I put in from AT2018zr
    """

    test_json = {
        "schema_version": {"value": "0", "comment": "Copied from tde.space"},
        "distance": [
            {
                "value": "0.354",
                "reference": [
                    "2011Natur.476..425Z",
                    "2012ApJ...748...36B",
                    "2012MNRAS.421.1942W",
                    "2013ApJ...767..152Z",
                    "2016MNRAS.462L..66Y",
                    "2018ApJ...854...86E",
                ],
                "computed": False,
                "default": True,
                "distance_type": "redshift",
            },
            {
                "value": "0.3543",
                "reference": ["2017ApJ...838..149A", "2011Sci...333..203B"],
                "computed": False,
                "default": True,
                "distance_type": "redshift",
            },
            {
                "value": "1942",
                "reference": [
                    "2017ApJ...838..149A",
                    "2016A&A...594A..13P",
                    "2011Sci...333..203B",
                    "2017ApJ...835...64G",
                ],
                "computed": False,
                "default": True,
                "distance_type": "luminosity",
                "unit": "Mpc",
            },
            {
                "value": "1434",
                "reference": [
                    "2017ApJ...838..149A",
                    "2016A&A...594A..13P",
                    "2011Sci...333..203B",
                    "2017ApJ...835...64G",
                ],
                "computed": False,
                "default": True,
                "distance_type": "comoving",
                "unit": "Mpc",
            },
        ],
        "filter_alias": [
            {"filter_key": "4.9GHz", "freq_eff": 4.9, "freq_units": "GHz"},
            {
                "filter_key": "0.3 - 2.0",
                "wave_eff": 0.7293188143129428,
                "wave_min": 0.6199209921660013,
                "wave_max": 4.132806614440009,
                "wave_units": "nm",
            },
            {"filter_key": "g", "wave_eff": 1, "wave_units": "nm"},
            {"filter_key": "r", "wave_eff": 1, "wave_units": "nm"},
        ],
        "reference_alias": [
            {
                "name": "2011Natur.476..425Z",
                "human_readable_name": "Zauderer et al. (2011)",
            },
            {
                "name": "2012ApJ...748...36B",
                "human_readable_name": "Berger et al. (2012)",
            },
            {
                "name": "2012MNRAS.421.1942W",
                "human_readable_name": "Wiersema et al. (2012)",
            },
            {
                "name": "2013ApJ...767..152Z",
                "human_readable_name": "Zauderer et al. (2013)",
            },
            {
                "name": "2016MNRAS.462L..66Y",
                "human_readable_name": "Yang et al. (2016)",
            },
            {
                "name": "2018ApJ...854...86E",
                "human_readable_name": "Eftekhari et al. (2018)",
            },
            {
                "name": "2017ApJ...838..149A",
                "human_readable_name": "Auchettl, Guillochon, & Ramirez-Ruiz (2017)",
            },
            {
                "name": "2016A&A...594A..13P",
                "human_readable_name": "Planck Collaboration et al. (2016)",
            },
            {
                "name": "2011Sci...333..203B",
                "human_readable_name": "Bloom et al. (2011)",
            },
            {
                "name": "2011ApJ...737..103S",
                "human_readable_name": "Schlafly & Finkbeiner (2011)",
            },
            {
                "name": "2017ApJ...835...64G",
                "human_readable_name": "Guillochon et al. (2017)",
            },
        ],
        "date_reference": [
            {
                "value": 55648.54,
                "date_format": "mjd",
                "reference": [
                    "2011Natur.476..425Z",
                    "2012ApJ...748...36B",
                    "2012MNRAS.421.1942W",
                    "2013ApJ...767..152Z",
                    "2016MNRAS.462L..66Y",
                    "2018ApJ...854...86E",
                ],
                "computed": False,
                "date_type": "discovery",
            },
            {
                "value": "2011-03-29",
                "date_format": "iso",
                "reference": ["2017ApJ...838..149A", "2017ApJ...835...64G"],
                "computed": False,
                "date_type": "discovery",
                "measurement_type": "discovery",
                "default": True,
            },
        ],
        "name": {
            "default_name": "Sw J1644+57",
            "alias": [
                {"value": "Sw J1644+57", "reference": ["Swift"]},
                {"value": "GRB 110328A", "reference": [["2011Sci...333..203B"]]},
                {
                    "value": "Swift J164449.3+573451",
                    "reference": [["2011Sci...333..203B"]],
                },
                {
                    "value": "Swift J1644+57",
                    "reference": [["2017ApJ...838..149A", "2011Sci...333..203B"]],
                },
            ],
        },
        "photometry": [
            {
                "reference": [
                    "2011Natur.476..425Z",
                    "2012ApJ...748...36B",
                    "2012MNRAS.421.1942W",
                    "2013ApJ...767..152Z",
                    "2016MNRAS.462L..66Y",
                    "2018ApJ...854...86E",
                ],
                "raw": [
                    0.25,
                    0.34,
                    0.34,
                    0.61,
                    0.82,
                    1.48,
                    1.47,
                    1.8,
                    2.1,
                    4.62,
                    4.84,
                    5.86,
                    9.06,
                    9.1,
                    9.1,
                    11.71,
                    12.93,
                    12.83,
                    13.29,
                    12.43,
                    12.17,
                    12.05,
                    12.24,
                    11.12,
                    8.9,
                    8.24,
                    8.63,
                    6.23,
                    4.21,
                    3.52,
                    2.34,
                    1.47,
                ],
                "raw_units": [
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                    "mJy",
                ],
                "date": [
                    55652.41,
                    55653.3,
                    55653.54,
                    55654.33,
                    55655.32,
                    55656.31,
                    55658.33,
                    55663.520000000004,
                    55671.32,
                    55684.4,
                    55699.19,
                    55716.15,
                    55743.18,
                    55760.16,
                    55775.05,
                    55792.16,
                    55812.92,
                    55823.01,
                    55845.950000000004,
                    55861.86,
                    55893.770000000004,
                    55951.49,
                    56032.46,
                    56102.200000000004,
                    56230.85,
                    56293.54,
                    56299.64,
                    56436.14,
                    56680.54,
                    56753.54,
                    57021.54,
                    57542.54,
                ],
                "date_format": [
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                    "mjd",
                ],
                "filter_key": [
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                    "4.9GHz",
                ],
                "computed": [
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
                "obs_type": [
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                    "radio",
                ],
                "upperlimit": [
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
                "corr_k": False,
                "corr_s": False,
                "corr_av": False,
                "corr_host": False,
                "corr_hostav": False,
            },
            {
                "raw": ["0.000427131", "0.000399543"],
                "raw_units": ["ct", "ct"],
                "value": ["3.05E-15", "2.86E-15"],
                "value_units": ["ergs/cm^2/s", "ergs/cm^2/s"],
                "filter_key": ["0.3 - 2.0", "0.3 - 2.0"],
                "obs_type": ["xray", "xray"],
                "date": ["56257.4", "57070.2"],
                "date_format": "MJD",
                "upperlimit": [True, True],
                "telescope": "Chandra",
                "corr_k": None,
                "corr_av": None,
                "corr_host": None,
                "corr_hostav": None,
                "reference": ["2017ApJ...838..149A"],
                "corr_s": None,
            },
            {
                "raw": ["0.0308041", "0.0993732", "0.0450751", "0.043094"],
                "raw_units": ["ct", "ct", "ct", "ct"],
                "value": ["8.75E-14", "1.66E-13", "7.53E-14", "7.20E-14"],
                "value_units": [
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                ],
                "filter_key": ["0.3 - 2.0", "0.3 - 2.0", "0.3 - 2.0", "0.3 - 2.0"],
                "obs_type": ["xray", "xray", "xray", "xray"],
                "date": ["48000", "48094.1", "48370", "48715.6"],
                "date_format": "MJD",
                "upperlimit": [True, True, True, True],
                "telescope": "ROSAT",
                "corr_k": None,
                "corr_av": None,
                "corr_host": None,
                "corr_hostav": None,
                "reference": ["2017ApJ...838..149A"],
                "corr_s": None,
            },
            {
                "raw": [
                    "1.0482",
                    "0.6817",
                    "0.101597",
                    "0.0362351",
                    "0.0225792",
                    "0.0122362",
                    "0.00300771",
                    "0.00126983",
                    "0.00156039",
                    "0.00173297",
                    "1.95E-03",
                    "1.86E-03",
                    "2.09E-03",
                    "1.97E-03",
                    "0.0019953",
                    "0.00171713",
                    "0.00184643",
                    "0.00252298",
                    "0.00236432",
                ],
                "raw_units": [
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                ],
                "value": [
                    "8.74E-12",
                    "5.53E-12",
                    "8.91E-13",
                    "3.04E-13",
                    "2.27E-13",
                    "9.46E-14",
                    "2.26E-14",
                    "2.21E-14",
                    "2.71E-14",
                    "3.00E-14",
                    "3.39E-14",
                    "3.23E-14",
                    "3.63E-14",
                    "3.42E-14",
                    "3.47E-14",
                    "2.99E-14",
                    "3.21E-14",
                    "4.38E-14",
                    "4.10E-14",
                ],
                "value_units": [
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                ],
                "filter_key": [
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                ],
                "obs_type": [
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                ],
                "date": [
                    "55649.6",
                    "55675.6",
                    "55749.3",
                    "55850.3",
                    "55947.6",
                    "56048.1",
                    "56153.6",
                    "56242.9",
                    "56335.7",
                    "56457",
                    "56551.70",
                    "56649.30",
                    "56751.40",
                    "56848.30",
                    "56945.1",
                    "57049",
                    "57151.3",
                    "57253.9",
                    "57318.4",
                ],
                "date_format": "MJD",
                "upperlimit": [
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
                "telescope": "Swift",
                "corr_k": None,
                "corr_av": None,
                "corr_host": None,
                "corr_hostav": None,
                "reference": ["2017ApJ...838..149A"],
                "corr_s": None,
            },
            {
                "raw": [
                    "0.0319464",
                    "11.3004",
                    "7.93976",
                    "0.928866",
                    "17.5",
                    "0.241757",
                    "1.72019",
                    "2.83",
                    "1.45692",
                    "1.51279",
                    "0.348197",
                    "0.614621",
                    "0.860492",
                    "0.0138642",
                    "1.08E-02",
                    "1.8",
                    "5.35E-01",
                ],
                "raw_units": [
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                    "ct",
                ],
                "value": [
                    "1.89E-14",
                    "9.03E-12",
                    "5.14E-12",
                    "6.93E-13",
                    "3.01E-11",
                    "1.53E-13",
                    "7.87E-13",
                    "4.86E-12",
                    "1.19E-13",
                    "4.07E-13",
                    "1.46E-13",
                    "3.64E-13",
                    "5.08E-13",
                    "8.22E-15",
                    "6.38E-15",
                    "3.09E-12",
                    "9.19E-13",
                ],
                "value_units": [
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                    "ergs/cm^2/s",
                ],
                "filter_key": [
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                    "0.3 - 2.0",
                ],
                "obs_type": [
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                    "xray",
                ],
                "date": [
                    "55651.4",
                    "55667.4",
                    "55681.3",
                    "55697.3",
                    "55697.6",
                    "55745.3",
                    "55757.1",
                    "55757.4",
                    "55769.1",
                    "55787.1",
                    "55801",
                    "55811",
                    "55836.9",
                    "56197.8",
                    "56205.8",
                    "56296.3",
                    "57231.4",
                ],
                "date_format": "MJD",
                "upperlimit": [
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
                "telescope": "XMM",
                "corr_k": None,
                "corr_av": None,
                "corr_host": None,
                "corr_hostav": None,
                "reference": ["2017ApJ...838..149A"],
                "corr_s": None,
            },
            {
                "raw": [
                    "21.44",
                    "21.13",
                    "21.65",
                    "21.06",
                    "20.85",
                    "20.47",
                    "18.85",
                    "18.84",
                    "18.77",
                    "18.64",
                    "18.64",
                    "18.16",
                    "17.98",
                    "17.59",
                    "17.59",
                ],
                "raw_units": "mag(AB)",
                "filter_key": [
                    "r",
                    "r",
                    "r",
                    "r",
                    "r",
                    "r",
                    "g",
                    "g",
                    "g",
                    "g",
                    "g",
                    "r",
                    "r",
                    "r",
                    "r",
                ],
                "obs_type": [
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                    "uvoir",
                ],
                "date": [
                    "58101.360",
                    "58105.310",
                    "58154.220",
                    "58156.230",
                    "58158.230",
                    "58160.230",
                    "58166.260",
                    "58167.180",
                    "58168.164",
                    "58168.176",
                    "58170.392",
                    "58182.190",
                    "58183.180",
                    "58190.154",
                    "58190.155",
                ],
                "date_format": "MJD",
                "upperlimit": [
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                ],
                "telescope": "P48",
                "corr_k": None,
                "corr_av": None,
                "corr_host": None,
                "corr_hostav": None,
                "reference": ["2018arXiv180902608V"],
                "corr_s": None,
            },
        ],
        "coordinate": [
            {
                "ra": "16 44 49.93130",
                "dec": "+57 34 59.6893",
                "epoch": "J2000",
                "system": "ICRS",
                "ra_units": "hourangle",
                "dec_units": "deg",
                "reference": [
                    "2011Natur.476..425Z",
                    "2012ApJ...748...36B",
                    "2012MNRAS.421.1942W",
                    "2013ApJ...767..152Z",
                    "2016MNRAS.462L..66Y",
                    "2018ApJ...854...86E",
                ],
                "computed": False,
                "default": True,
                "uuid": "b98cdbe9-fff7-415b-a6d6-78cbd88544b8",
                "coordinate_type": "equatorial",
            },
            {
                "reference": "b98cdbe9-fff7-415b-a6d6-78cbd88544b8",
                "computed": True,
                "coordinate_type": "galactic",
                "l": 86.71123273634593,
                "b": 39.44122227143398,
                "l_units": "deg",
                "b_units": "deg",
            },
            {
                "ra": "16:44:49.92",
                "dec": "+57:34:58.8",
                "ra_units": "hour",
                "dec_units": "deg",
                "coordinate_type": "equitorial",
            },
        ],
        "classification": [
            {
                "object_class": "TDE",
                "confidence": 1.0,
                "reference": [
                    "2011Natur.476..425Z",
                    "2011Sci...333..203B",
                    "2012ApJ...748...36B",
                    "2012MNRAS.421.1942W",
                    "2013ApJ...767..152Z",
                    "2016MNRAS.462L..66Y",
                    "2018ApJ...854...86E",
                ],
                "default": True,
            }
        ],
        "host": [
            {
                "host_name": "Swift J164449.3+573451",
                "host_ra": 251.20416666666665,
                "host_dec": 57.58083333333334,
                "host_ra_units": "deg",
                "host_dec_units": "deg",
                "reference": [
                    "2023PASP..135c4101G",
                    "2011Sci...333..199L",
                    "2011Sci...333..203B",
                    "2017ApJ...838..149A",
                ],
            }
        ],
    }

    return test_json
