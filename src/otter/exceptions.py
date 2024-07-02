"""
Custom exceptions for otter
"""

from __future__ import annotations


class FailedQueryError(ValueError):
    """
    Exception thrown when the users query does not return any results.
    """

    def __str__(self):
        txt = "You're query/search did not return any results! "
        txt += "Try again with different parameters!"
        return txt


class IOError(ValueError):
    """
    Exception thrown when the input or output argument/value is not the correct type.
    """

    pass


class OtterLimitationError(Exception):
    """
    Exception thrown when the user requests something that is currently not supported
    by the API.
    """

    def __init__(self, msg):
        self.msg = "Current Limitation Found: " + msg

    def __str__(self):
        return self.msg


class TransientMergeError(Exception):
    """
    Exception thrown when the Transient objects can not be combined as expected.
    """

    pass


class OtterNotImplementedError(Exception):
    """
    Exception thrown by methods that are not yet implemented
    """

    def __str__(self):
        return "Not yet implemented! Please consider contributing on our GitHub!"


class MissingEnvVarError(Exception):
    """
    Thrown when the querying code expects an Environment variable and doesn't get one
    """

    def __init__(self, v, url=None):
        self.v = v
        self.url = url

    def __str__(self):
        msg = "You are missing the following environment variables:\n"
        for vv in self.v:
            msg += f"\t- {vv}\n"
        msg += "Please make an account and/or set these variables!\n"
        if self.url is not None:
            msg += f"See: {self.url}"

        return msg
