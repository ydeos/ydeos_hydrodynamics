#!/usr/bin/env python
# coding: utf-8

r"""Tests for the mit.py module."""

import pytest

from ydeos_hydrodynamics.mit import hull_residuary_resistance_mit


def test_happy_path():
    r"""Test a normal call."""
    lwl = 5.
    Tc = 0.25
    force = hull_residuary_resistance_mit(boatspeed=2.,
                                          Vc=1.,
                                          lwl=lwl,
                                          bwl=0.8,
                                          Tc=Tc,
                                          coe=None)
    assert force.fx < 0
    assert force.fy == 0
    assert force.fz == 0
    assert force.px == lwl / 2
    assert force.py == 0
    assert force.pz == -Tc / 3


def test_wrong_input():
    r"""Test bad input."""
    with pytest.raises(ValueError):
        hull_residuary_resistance_mit(boatspeed=2.,
                                      Vc=0.,
                                      lwl=5.,
                                      bwl=0.8,
                                      Tc=0.25,
                                      coe=None)
    with pytest.raises(ValueError):
        hull_residuary_resistance_mit(boatspeed=2.,
                                      Vc=1.,
                                      lwl=0.,
                                      bwl=0.8,
                                      Tc=0.25,
                                      coe=None)
    with pytest.raises(ValueError):
        hull_residuary_resistance_mit(boatspeed=2.,
                                      Vc=1.,
                                      lwl=5.,
                                      bwl=0.,
                                      Tc=0.25,
                                      coe=None)
    with pytest.raises(ValueError):
        hull_residuary_resistance_mit(boatspeed=2.,
                                      Vc=1.,
                                      lwl=5.,
                                      bwl=0.8,
                                      Tc=0.,
                                      coe=None)
    with pytest.raises(ValueError):
        hull_residuary_resistance_mit(boatspeed=2.,
                                      Vc=1.,
                                      lwl=5.,
                                      bwl=0.8,
                                      Tc=0.25,
                                      coe=None,
                                      rho_water=900.)
