"""Tests for the ToroFit package."""

import torofit
from torofit import version_info
from torofit.core import version_info as core_version_info


def test_version():
    assert torofit.__version__ == "0.1.0"


def test_version_info():
    info = version_info()
    assert info["name"] == "torofit"
    assert info["version"] == "0.1.0"


def test_version_info_from_core():
    info = core_version_info()
    assert info["name"] == "torofit"
    assert info["version"] == "0.1.0"
