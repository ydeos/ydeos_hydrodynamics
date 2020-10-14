# coding: utf-8

r"""Tests for the domain.py module"""

from ydeos_hydrodynamics.domain import check_domain


def test_domain_inside():
    r"""The value is inside of the domain."""
    assert check_domain("Whatever", "value", lower=1, upper=2, value=3) is False


def test_domain_outside():
    r"""The value is outside of the domain."""
    assert check_domain("Whatever", "value", lower=1, upper=3, value=2) is True
