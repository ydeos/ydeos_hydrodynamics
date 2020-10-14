#!/usr/bin/env python
# coding: utf-8

r"""Tests for the viscous.py module."""

import pytest
from ydeos_hydrodynamics.viscous import appendage_viscous, ballast_viscous, \
    hull_viscous


avra = appendage_viscous(boatspeed=3.5,
                         heel_angle=0.,
                         coe=(0.5, 0., -0.2),
                         average_chord=0.08,
                         thickness_to_chord=0.06,
                         wetted_area=2 * 0.02)


def test_known_value_appendage_viscous():
    r"""p 64 of Principles of Yacht Design - Larsson / Eliasson.

    thickness_to_chord=0.,  # form factor is 1 in book

    """
    avra = appendage_viscous(boatspeed=3.5,
                             heel_angle=0.,
                             coe=(0.5, 0., -0.2),
                             average_chord=1.45,
                             thickness_to_chord=0.,
                             wetted_area=4.4)
    # thickness_to_chord=0., # form factor is 1 in book
    avra_2e5 = appendage_viscous(boatspeed=3.5,
                                 heel_angle=0.,
                                 coe=(0.5, 0., -0.2),
                                 average_chord=1.45,
                                 thickness_to_chord=0.,
                                 wetted_area=4.4,
                                 transition_reynolds_number=2e5)
    # ittc57
    # assert abs(avra.force(boatspeed = 3.5).force[0]), 93.64061076321511)
    # prandtl
    assert abs(avra.fx) == 83.79539691196497
    # avra_2e5 has more turbulent flow than avra
    assert abs(avra.fx) < abs(avra_2e5.fx)


def test_negative_avg_chord_appendage_viscous():
    r"""Negative average chord."""
    with pytest.raises(ValueError):
        appendage_viscous(boatspeed=3.5,
                          heel_angle=0.,
                          coe=(0.5, 0., -0.2),
                          average_chord=-0.08,
                          thickness_to_chord=0.06,
                          wetted_area=0.02)


def test_coe_shift_appendage_viscous():
    r"""Test position change of the centre of effort."""
    force = appendage_viscous(boatspeed=1.,
                              heel_angle=30.,
                              coe=(0.5, 0., -0.2),
                              average_chord=0.08,
                              thickness_to_chord=0.06,
                              wetted_area=2 * 0.02)
    assert -0.1 - 1e-6 <= force.py <= -0.1 + 1e-6
    force = appendage_viscous(boatspeed=1.,
                              heel_angle=-30.,
                              coe=(0.5, 0., -0.2),
                              average_chord=0.08,
                              thickness_to_chord=0.06,
                              wetted_area=2 * 0.02)
    assert 0.1 - 1e-6 <= force.py <= 0.1 + 1e-6

    force = appendage_viscous(boatspeed=1.,
                              heel_angle=60.,
                              coe=(0.5, 0., -0.2),
                              average_chord=0.08,
                              thickness_to_chord=0.06,
                              wetted_area=2 * 0.02)
    assert -0.1 - 1e-6 <= force.pz <= -0.1 + 1e-6
    force = appendage_viscous(boatspeed=1.,
                              heel_angle=-60.,
                              coe=(0.5, 0., -0.2),
                              average_chord=0.08,
                              thickness_to_chord=0.06,
                              wetted_area=2 * 0.02)
    assert -0.1 - 1e-6 <= force.pz <= -0.1 + 1e-6


def test_resistance_sign_and_symmetry_appendage_viscous():
    r"""Test signs and symmetry."""
    force = appendage_viscous(boatspeed=1.,
                              heel_angle=30.,
                              coe=(0.5, 0., -0.2),
                              average_chord=0.08,
                              thickness_to_chord=0.06,
                              wetted_area=2 * 0.02)
    force1 = appendage_viscous(boatspeed=-1.,
                               heel_angle=30.,
                               coe=(0.5, 0., -0.2),
                               average_chord=0.08,
                               thickness_to_chord=0.06,
                               wetted_area=2 * 0.02)
    assert force.fx == -force1.fx
    assert force.fx < 0


bvra = ballast_viscous(boatspeed=1.,
                       heel_angle=0.,
                       coe=(0.5, 0., -0.2),
                       bulb_length=0.37,
                       bulb_diameter=0.037,
                       bulb_wetted_area=0.02)


def test_known_value_ballast_viscous():
    r"""Tests with different transition Reynolds."""
    bvra_2e5 = ballast_viscous(boatspeed=1.,
                               heel_angle=0.,
                               coe=(0.5, 0., -0.2),
                               bulb_length=0.37,
                               bulb_diameter=0.037,
                               bulb_wetted_area=0.02,
                               transition_reynolds_number=2e5)
    assert abs(bvra_2e5.fx) > abs(bvra.fx)


def test_negative_length_ballast_viscous():
    r"""Negative length."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=-0.37,
                        bulb_diameter=0.037,
                        bulb_wetted_area=0.02)


def test_coe_shift_ballast_viscous():
    r"""Test position change of centre of effort."""
    force = ballast_viscous(boatspeed=1.,
                            heel_angle=30.,
                            coe=(0.5, 0., -0.2),
                            bulb_length=0.37,
                            bulb_diameter=0.037,
                            bulb_wetted_area=0.02)
    assert -0.1 - 1e-6 <= force.py <= -0.1 + 1e-6
    force = ballast_viscous(boatspeed=1.,
                            heel_angle=-30.,
                            coe=(0.5, 0., -0.2),
                            bulb_length=0.37,
                            bulb_diameter=0.037,
                            bulb_wetted_area=0.02)
    assert 0.1 - 1e-6 <= force.py <= 0.1 + 1e-6

    force = ballast_viscous(boatspeed=1.,
                            heel_angle=60.,
                            coe=(0.5, 0., -0.2),
                            bulb_length=0.37,
                            bulb_diameter=0.037,
                            bulb_wetted_area=0.02)
    assert -0.1 - 1e-6 <= force.pz <= -0.1 + 1e-6
    force = ballast_viscous(boatspeed=1.,
                            heel_angle=-60.,
                            coe=(0.5, 0., -0.2),
                            bulb_length=0.37,
                            bulb_diameter=0.037,
                            bulb_wetted_area=0.02)
    assert -0.1 - 1e-6 <= force.pz <= -0.1 + 1e-6


def test_resistance_sign_and_symmetry_ballast_viscous():
    r"""Test sign and symmetry."""
    force = ballast_viscous(boatspeed=1.,
                            heel_angle=30.,
                            coe=(0.5, 0., -0.2),
                            bulb_length=0.37,
                            bulb_diameter=0.037,
                            bulb_wetted_area=0.02)
    force1 = ballast_viscous(boatspeed=-1.,
                             heel_angle=30.,
                             coe=(0.5, 0., -0.2),
                             bulb_length=0.37,
                             bulb_diameter=0.037,
                             bulb_wetted_area=0.02)
    assert force.fx == -force1.fx
    assert force.fx < 0


hfra = hull_viscous(boatspeed=3.5,
                    coe=(0.5, 0., -0.02),
                    lwl=10., Sc=25.2,
                    rho_water=1025.)


def test_known_value_hull_viscous():
    r"""  p 64 of Principles of Yacht Design - Larsson / Eliasson.

    408.8 found instead of 408.6 in book

    """
    # with ittc57
    # assert abs(hfra.fx) == 408.81620236106744
    # with prandtl method"""
    assert abs(hfra.fx) == 397.13619865572065


def test_negative_lwl_hull_viscous():
    r"""Test a negative length at waterline."""
    # Negative lwl
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=-10.,
                     Sc=25.2,
                     rho_water=1025.)


def test_resistance_sign_and_symmetry_hull_viscous():
    r"""Test the resistance sign and symmetry."""
    force = hull_viscous(boatspeed=1.,
                         coe=(0.5, 0., -0.02),
                         lwl=10.,
                         Sc=25.2,
                         rho_water=1025.)
    force_backwards = hull_viscous(boatspeed=-1.,
                                   coe=(0.5, 0., -0.02),
                                   lwl=10.,
                                   Sc=25.2,
                                   rho_water=1025.)
    assert force.fx == -force_backwards.fx
    assert force.fx < 0
