"""
Custom exceptions for otter
"""


class FailedQueryError(ValueError):
    def __str__(self):
        txt = "You're query/search did not return any results! "
        txt += "Try again with different parameters!"
        return txt


class IOError(ValueError):
    pass


class OtterLimitationError(Exception):
    def __init__(self, msg):
        self.msg = "Current Limitation Found: " + msg

    def __str__(self):
        return self.msg


class TransientMergeError(Exception):
    pass
