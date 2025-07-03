"""
Test the utility functions used throughout the OTTER API
"""

import json
from copy import deepcopy
from otter import util
from astropy import units as u
import pytest


def test_filter_to_obstype():
    """
    Test the converion from a filtername to a obstype variable. This is used a lot
    during the data cleaning process
    """

    assert util.filter_to_obstype("FUV") == "uvoir"
    assert util.filter_to_obstype("r") == "uvoir"

    with pytest.raises(Exception):
        util.filter_to_obstype("some fake filter name")


def test_freq_to_obstype():
    """
    Test the frequency conversion to an obstype
    """
    assert util.freq_to_obstype(1 * u.GHz) == "radio"


def test_wave_to_obstype():
    assert util.wave_to_obstype(100 * u.nm) == "uvoir"
    assert util.wave_to_obstype(1 * u.m) == "radio"
    assert util.wave_to_obstype(1e-3 * u.nm) == "xray"


def test_clean_schema():
    """
    Also used a lot during the data cleaning process. This function tests the
    clean_schema method.
    """

    subschema = deepcopy(util.distance_schema)

    subschema["value"] = "foo"
    subschema["reference"] = "bar"

    cleaned_schema = util.clean_schema(subschema)
    assert "unit" not in cleaned_schema
    assert "value" in cleaned_schema
    assert "reference" in cleaned_schema
    assert "error" not in cleaned_schema
    assert cleaned_schema["value"] == "foo"
    assert cleaned_schema["reference"] == "bar"


def test_bibcode_to_hrn():
    """
    Test the bibcode to hrn that only relies on a local file
    """
    bibcode_to_hrn = {"foo": "bar", "xyz": "wqr"}

    with open("reference_map_local.json", "w") as f:
        json.dump(bibcode_to_hrn, f, indent=4)

    # now we can actually test the local version
    bibcodes, hrns = util.bibcode_to_hrn("foo")
    assert bibcodes == ["foo"]
    assert hrns == ["bar"]

    bibcodes, hrns = util.bibcode_to_hrn(["foo", "xyz"])
    assert bibcodes == ["foo", "xyz"]
    assert hrns == ["bar", "wqr"]


def test_freq_to_band():
    """
    Test the conversion from a frequency quantity to a single band
    """

    res = util.freq_to_band(1.4 * u.GHz)
    assert res == "L"

    res = util.freq_to_band(10 * u.GHz)
    assert res == "X"

    with pytest.raises(ValueError):
        res = util.freq_to_band(1 * u.Hz)


def test_freqlist_to_band():
    res = util.freqlist_to_band([1, 10], ["GHz", "GHz"])
    assert res[0] == "UHF" and res[1] == "X"
