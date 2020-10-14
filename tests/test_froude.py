#!/usr/bin/env python
# coding: utf-8

r"""Tests for the froude.py module"""

import pytest
from ydeos_hydrodynamics.froude import froude_number


def test_known_value_froude_number():
    r"""Test against a known value."""
    target = 0.15
    tolerance = 0.001
    assert target - tolerance <= froude_number(speed=0.594, lwl=1.6) <= target + tolerance


def test_symmetry_froude_number():
    r"""Froude number symmetry forward and backward.

    Test the froude number is the same for positive and negative speeds with
    the same absolute value

    """
    assert froude_number(speed=0.594, lwl=1.6) == froude_number(speed=-0.594, lwl=1.6)


def test_wrong_input_froude_number():
    r"""Test wrong parameters."""
    with pytest.raises(ZeroDivisionError):
        _ = froude_number(speed=1, lwl=0)
