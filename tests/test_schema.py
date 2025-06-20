"""
Test the schema definitions and validators
"""

from otter.schema import (
    CoordinateSchema,
    DistanceSchema,
    PhotometrySchema,
    HostSchema,
    OtterSchema,
)
from pydantic import ValidationError
import pytest


def test_coordinate_validator():
    # first, some that should succeed
    case1 = dict(reference="foo", ra=1, dec=1, ra_units="deg", dec_units="deg")
    cs1 = CoordinateSchema(**case1)
    assert cs1.ra == 1 and cs1.dec == 1

    case2 = dict(reference="foo", l=1, b=1, l_units="deg", b_units="deg")
    cs1 = CoordinateSchema(**case2)
    assert cs1.l == 1 and cs1.b == 1

    case3 = dict(reference="foo", lon=1, lat=1, lon_units="deg", lat_units="deg")
    cs1 = CoordinateSchema(**case3)
    assert cs1.lon == 1 and cs1.lat == 1

    # now some cases that won't work
    del case1["dec_units"]
    with pytest.raises(ValueError):
        CoordinateSchema(**case1)

    del case1["ra_units"]
    with pytest.raises(ValueError):
        CoordinateSchema(**case1)

    del case2["b_units"]
    with pytest.raises(ValueError):
        CoordinateSchema(**case2)

    del case2["l_units"]
    with pytest.raises(ValueError):
        CoordinateSchema(**case2)

    del case3["lat_units"]
    with pytest.raises(ValueError):
        CoordinateSchema(**case3)

    del case3["lon_units"]
    with pytest.raises(ValueError):
        CoordinateSchema(**case3)

    case4 = {"reference": "foo"}
    with pytest.raises(ValueError):
        CoordinateSchema(**case4)


def test_distance_validator():
    case1 = dict(value=1, reference="foo", distance_type="redshift")

    case2 = dict(value=1, reference="foo", distance_type="luminosity", unit="cm")

    cs1 = DistanceSchema(**case1)
    assert cs1.value == 1
    assert cs1.reference == "foo"
    assert cs1.unit is None

    cs2 = DistanceSchema(**case2)
    assert cs2.value == 1
    assert cs2.reference == "foo"
    assert cs2.distance_type == "luminosity"
    assert cs2.unit == "cm"

    del case2["unit"]
    with pytest.raises(ValueError):
        DistanceSchema(**case2)


def test_photometry_validators():
    case1 = dict(
        reference="foo",
        raw=[1, 2],
        raw_err=[3, 4],
        raw_units="mJy",
        filter_key="X",
        obs_type="radio",
        date=[1, 2],
        date_format="mjd",
    )

    PhotometrySchema(**case1)

    case2 = dict(
        reference="foo",
        raw=1,
        raw_err=3,
        raw_units="mJy",
        filter_key="X",
        obs_type="radio",
        date=1,
        date_format="mjd",
    )

    with pytest.raises(ValueError):
        PhotometrySchema(**case2)


def test_host_validator():
    case1 = dict(reference="foo", ra=1, dec=1, ra_units="deg")

    with pytest.raises(ValueError):
        HostSchema(**case1)

    del case1["ra_units"]
    with pytest.raises(ValueError):
        HostSchema(**case1)

    case2 = dict(reference="foo", host_ra_units="deg", host_dec_units="deg")

    with pytest.raises(ValueError):
        HostSchema(**case2)

    case3 = dict(
        reference="foo",
        host_ra=1,
        host_ra_units="deg",
        host_dec_units="deg",
        host_name="foo",
    )

    with pytest.raises(ValueError):
        HostSchema(**case3)

    case4 = dict(
        reference="foo",
        host_ra=1,
        host_dec=1,
        host_ra_units="deg",
        host_dec_units="deg",
        host_name="foo",
    )

    case5 = dict(reference="foo", host_name="foo")

    case6 = dict(
        host_ra=1,
        host_dec=1,
        host_ra_units="deg",
        host_dec_units="deg",
    )

    HostSchema(**case4)
    HostSchema(**case5)

    with pytest.raises(ValidationError):
        HostSchema(**case6)


def test_otter_validator():
    case1 = dict(
        name=dict(default_name="foo", alias=[dict(value="foo", reference="foo")]),
        coordinate=[
            dict(ra=1, dec=1, ra_units="deg", dec_units="deg", reference="foo")
        ],
        photometry=[
            dict(
                reference="foo",
                raw=[1, 2],
                raw_err=[3, 4],
                raw_units="mJy",
                filter_key="X",
                obs_type="radio",
                date=[1, 2],
                date_format="mjd",
            )
        ],
        reference_alias=[{"name": "foo", "human_readable_name": "foo"}],
    )

    with pytest.raises(ValueError):
        OtterSchema(**case1)

    case1["filter_alias"] = [
        dict(filter_key="X", filter_name="X", wave_eff=1, wave_units="nm")
    ]

    OtterSchema(**case1)
