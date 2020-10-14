#!/usr/bin/env python
# coding: utf-8

r"""Tests for the righting_moment.py module."""

import pytest
from ydeos_hydrodynamics.righting_moment import rm_estimate_gerritsma, \
    trm_estimate
from ydeos_forces.forces import Force, SystemOfForces

displacement = 4.0
lwl = 1.0
bwl = 0.2
Tc = 0.05
Gdwl = -0.2
rho_water = 1025.

rm = rm_estimate_gerritsma(boatspeed=1.,
                           heel_angle=30.,
                           displacement=displacement,
                           lwl=lwl,
                           bwl=bwl,
                           Tc=Tc,
                           Gdwl=Gdwl,
                           rho_water=rho_water)


def test_known_value_rm_estimate_gerritsma():
    r"""No known value."""
    assert True


def test_wrong_input_rm_estimate_gerritsma():
    r"""Problematic input for heeling righting moment."""
    # negative lwl
    with pytest.raises(ValueError):
        rm_estimate_gerritsma(boatspeed=1.,
                              heel_angle=30.,
                              displacement=displacement,
                              lwl=-lwl,
                              bwl=bwl,
                              Tc=Tc,
                              Gdwl=Gdwl,
                              rho_water=rho_water)
    # zero displacement
    with pytest.raises(ValueError):
        rm_estimate_gerritsma(boatspeed=1.,
                              heel_angle=30.,
                              displacement=0.,
                              lwl=lwl,
                              bwl=bwl,
                              Tc=Tc,
                              Gdwl=Gdwl,
                              rho_water=rho_water)
    # zero bwl
    with pytest.raises(ValueError):
        rm_estimate_gerritsma(boatspeed=1.,
                              heel_angle=30.,
                              displacement=displacement,
                              lwl=lwl,
                              bwl=0.,
                              Tc=Tc,
                              Gdwl=Gdwl,
                              rho_water=rho_water)
    # zero Tc
    with pytest.raises(ValueError):
        rm_estimate_gerritsma(boatspeed=1.,
                              heel_angle=30.,
                              displacement=displacement,
                              lwl=lwl,
                              bwl=bwl,
                              Tc=0.,
                              Gdwl=Gdwl,
                              rho_water=rho_water)
    # problematic rho_water
    with pytest.raises(ValueError):
        rm_estimate_gerritsma(boatspeed=1.,
                              heel_angle=30.,
                              displacement=displacement,
                              lwl=lwl,
                              bwl=bwl,
                              Tc=Tc,
                              Gdwl=Gdwl,
                              rho_water=930.)


def test_rm_sign_and_symmetry_rm_estimate_gerritsma():
    r"""Test righting moment sign."""
    s = SystemOfForces()
    for f in rm:
        s.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))
    assert s.mx > 0
    assert s.my == 0.
    assert s.mz == 0.


def test_rm_negative_heel_rm_estimate_gerritsma():
    r"""Test righting moment with negative heel."""
    s1 = SystemOfForces()
    for f in rm_estimate_gerritsma(boatspeed=1.,
                                   heel_angle=-30.,
                                   displacement=displacement,
                                   lwl=lwl,
                                   bwl=bwl,
                                   Tc=Tc,
                                   Gdwl=Gdwl,
                                   rho_water=rho_water):
        s1.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))
    assert s1.mx < 0
    assert s1.my == 0.
    assert s1.mz == 0.


def test_rm_going_backwards_rm_estimate_gerritsma():
    r"""Test righting moment going backwards."""
    s = SystemOfForces()
    for f in rm_estimate_gerritsma(boatspeed=-1.,
                                   heel_angle=30.,
                                   displacement=displacement,
                                   lwl=lwl,
                                   bwl=bwl,
                                   Tc=Tc,
                                   Gdwl=Gdwl,
                                   rho_water=rho_water):
        s.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))
    assert s.mx > 0
    assert s.my == 0.
    assert s.mz == 0.


def test_rm_zero_heel_rm_estimate_gerritsma():
    r"""Test righting moment and sum of forces are 0 when upright."""
    s = SystemOfForces()
    for f in rm_estimate_gerritsma(boatspeed=1.,
                                   heel_angle=0.,
                                   displacement=displacement,
                                   lwl=lwl,
                                   bwl=bwl,
                                   Tc=Tc,
                                   Gdwl=Gdwl,
                                   rho_water=rho_water):
        s.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))
    assert s.mx == 0.
    assert s.my == 0.
    assert s.mz == 0.

    assert s.x == 0.
    assert s.y == 0.
    assert s.z == 0.


displacement = 4.0
lwl = 1.0
bwl = 0.2
Tc = 0.05
Gdwl = -0.2
rho_water = 1025.

trm = trm_estimate(trim_angle=3.,
                   displacement=displacement,
                   lwl=lwl,
                   bwl=bwl,
                   Tc=Tc,
                   Gdwl=Gdwl,
                   rho_water=rho_water)


def test_known_value_trm_estimate():
    r"""No known value."""
    pass


def test_negative_lwl_trm_estimate():
    r"""Negative lwl."""
    with pytest.raises(ValueError):
        trm_estimate(trim_angle=3.,
                     displacement=displacement,
                     lwl=-lwl,
                     bwl=bwl,
                     Tc=Tc,
                     Gdwl=Gdwl,
                     rho_water=rho_water)


def test_trm_sign_and_symmetry_trm_estimate():
    r"""Test trim righting moment sign."""
    s = SystemOfForces()
    for f in trm:
        s.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))
    assert s.my > 0
    assert s.mx == 0.
    assert s.mz == 0.


def test_trm_negative_pitch_trm_estimate():
    r"""Test moments with negative trim."""
    s1 = SystemOfForces()
    for f in trm_estimate(trim_angle=-3.,
                          displacement=displacement,
                          lwl=lwl,
                          bwl=bwl,
                          Tc=Tc,
                          Gdwl=Gdwl,
                          rho_water=rho_water):
        s1.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))
    assert s1.mx == 0.
    assert s1.my < 0
    assert s1.mz == 0.


# def test_trm_going_backwards():
#     s = SystemOfForces()
#     for f in trm.trim_righting_moment(boatspeed=-1., trim_angle=3.):
#         s.add_force(f)
#     assert s.moment[1] > 0
#     assert s.moment[0] == 0.
#     assert s.moment[2] == 0.


def test_trm_zero_heel_moving_trm_estimate():
    r"""Check moments and sum pf forces is 0 when upright."""
    s = SystemOfForces()
    for f in trm_estimate(trim_angle=0.,
                          displacement=displacement,
                          lwl=lwl,
                          bwl=bwl,
                          Tc=Tc,
                          Gdwl=Gdwl):
        s.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))
    assert s.mx == 0.
    assert s.my == 0.
    assert s.mz == 0.

    assert s.x == 0.
    assert s.y == 0.
    assert s.z == 0.


def test_wrong_input_trm_estimate():
    r"""Problematic input for trim righting moment."""
    # negative lwl
    with pytest.raises(ValueError):
        trm_estimate(trim_angle=30.,
                     displacement=displacement,
                     lwl=-lwl,
                     bwl=bwl,
                     Tc=Tc,
                     Gdwl=Gdwl,
                     rho_water=rho_water)
    # zero displacement
    with pytest.raises(ValueError):
        trm_estimate(trim_angle=30.,
                     displacement=0.,
                     lwl=lwl,
                     bwl=bwl,
                     Tc=Tc,
                     Gdwl=Gdwl,
                     rho_water=rho_water)
    # zero bwl
    with pytest.raises(ValueError):
        trm_estimate(trim_angle=30.,
                     displacement=displacement,
                     lwl=lwl,
                     bwl=0.,
                     Tc=Tc,
                     Gdwl=Gdwl,
                     rho_water=rho_water)
    # zero Tc
    with pytest.raises(ValueError):
        trm_estimate(trim_angle=30.,
                     displacement=displacement,
                     lwl=lwl,
                     bwl=bwl,
                     Tc=0.,
                     Gdwl=Gdwl,
                     rho_water=rho_water)
    # problematic rho_water
    with pytest.raises(ValueError):
        trm_estimate(trim_angle=30.,
                     displacement=displacement,
                     lwl=lwl,
                     bwl=bwl,
                     Tc=Tc,
                     Gdwl=Gdwl,
                     rho_water=930.)
