#!/usr/bin/env python
# coding: utf-8

r"""Tests for the geometry.py module"""

from ydeos_hydrodynamics.geometry import shift_coe_around_x


def test_shift_coe_around_x():
    r"""Test COE rotation around the x axis"""
    assert shift_coe_around_x((1., 0., 0.), heel_angle=0.) == (1., 0., 0.)
    assert shift_coe_around_x((1., 1., 1.), heel_angle=0.) == (1., 1., 1.)
    assert shift_coe_around_x((1., 0., 0.), heel_angle=90.) == (1., 0., 0.)
    assert shift_coe_around_x((1., 1., 0.), heel_angle=90.)[0] == 1.
    assert shift_coe_around_x((1., 1., 0.), heel_angle=90.)[1] <= 1e-10
    assert shift_coe_around_x((1., 1., 0.), heel_angle=90.)[2] == -1.
