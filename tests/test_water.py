#!/usr/bin/env python
# coding: utf-8

r"""Tests for the water.py module."""

from ydeos_hydrodynamics.water import water_fresh_density, \
    water_fresh_dynamic_viscosity, water_fresh_kinematic_viscosity, \
    water_sea_density, water_sea_dynamic_viscosity, \
    water_sea_kinematic_viscosity


def test_water_fresh_density():
    r"""Test happy path."""
    assert water_fresh_density(10.) == 999.70


def test_water_sea_density():
    r"""Test happy path."""
    assert water_sea_density(10.) == 1027.


def test_water_fresh_dynamic_viscosity():
    r"""Test happy path."""
    assert water_fresh_dynamic_viscosity(10.) == 0.0013060


def test_water_sea_dynamic_viscosity():
    r"""Test happy path."""
    assert water_sea_dynamic_viscosity(10.) == 0.00141


def test_water_fresh_kinematic_viscosity():
    r"""Test happy path."""
    assert water_fresh_kinematic_viscosity(10.) == 1.3063919175752727e-06


def test_water_sea_kinematic_viscosity():
    r"""Test happy path."""
    assert water_sea_kinematic_viscosity(10.) == 1.3729308666017527e-06
