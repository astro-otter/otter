"""
Test the utility functions used throughout the OTTER API
"""

from copy import deepcopy
from otter import util


def test_filter_to_obstype():
    """
    Test the converion from a filtername to a obstype variable. This is used a lot
    during the data cleaning process
    """

    assert util.filter_to_obstype("FUV") == "uvoir"
    assert util.filter_to_obstype("r") == "uvoir"


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
    Since GitHub can't have an ADS API token we cant really implement these tests
    until we figure it out...
    """
    pass
