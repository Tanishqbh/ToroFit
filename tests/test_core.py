"""Tests for the ToroFit package."""

import torofit
from torofit import greet
from torofit.core import greet as core_greet


def test_version():
    assert torofit.__version__ == "0.1.0"


def test_greet_default():
    assert greet() == "Hello, World! Welcome to ToroFit."


def test_greet_with_name():
    assert greet("Athlete") == "Hello, Athlete! Welcome to ToroFit."


def test_greet_from_core():
    assert core_greet("Runner") == "Hello, Runner! Welcome to ToroFit."
