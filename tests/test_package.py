"""
Tests that the package actually imports and has a version
"""

import otter


def test_version():
    assert otter.__version__, "No version found for the otter package!"
