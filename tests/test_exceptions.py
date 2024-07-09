"""
Make sure all of the exceptions don't throw silly errors
"""

import pytest
from otter import exceptions as exc


def test_failed_query_error():
    with pytest.raises(exc.FailedQueryError):
        raise exc.FailedQueryError()

    txt = "You're query/search did not return any results! "
    txt += "Try again with different parameters!"

    assert str(exc.FailedQueryError()) == txt


def test_ioerror():
    with pytest.raises(exc.IOError):
        raise exc.IOError


def test_otter_limitation_error():
    with pytest.raises(exc.OtterLimitationError):
        raise exc.OtterLimitationError("")

    assert str(exc.OtterLimitationError("foo")) == "Current Limitation Found: foo"


def test_transient_merge_error():
    with pytest.raises(exc.TransientMergeError):
        raise exc.TransientMergeError


def test_otter_not_implemented_error():
    with pytest.raises(exc.OtterNotImplementedError):
        raise exc.OtterNotImplementedError

    assert str(exc.OtterNotImplementedError("foo")) == str(
        exc.OtterNotImplementedError()
    )


def test_missing_env_var_error():
    with pytest.raises(exc.MissingEnvVarError):
        raise exc.MissingEnvVarError(["val"])

    e = exc.MissingEnvVarError(["foo"])
    assert e.v == ["foo"]
    assert e.url is None
    assert str(e)[-1] == "\n"

    e = exc.MissingEnvVarError(["foo"], url="fakeurl.com")
    assert e.v == ["foo"]
    assert e.url == "fakeurl.com"
    assert str(e)[-1] == "m"
