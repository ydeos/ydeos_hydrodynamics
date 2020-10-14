#!/usr/bin/env python
# coding: utf-8

r"""Tests for the wetted_area.py module."""

from ydeos_hydrodynamics.wetted_area import wetted_area


def test_wetted_area_happy_path():
    tol = 1e-3
    expected = 20.2584
    wsa_heeled = wetted_area(20., bwl=3., Tc=1., Cm=0.75, heel=20.)
    assert expected - tol <= wsa_heeled <= expected + tol


def test_negative_heel():
    r"""Check the heel sign has no influence"""
    wsa_heel_plus = wetted_area(20., bwl=3., Tc=1., Cm=0.75, heel=20.)
    wsa_heel_minus = wetted_area(20., bwl=3., Tc=1., Cm=0.75, heel=-20.)
    assert wsa_heel_plus == wsa_heel_minus
