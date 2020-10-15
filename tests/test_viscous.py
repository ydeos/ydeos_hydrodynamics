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


# Wrong input hull


def test_negative_lwl_hull_viscous():
    r"""Test a negative length at waterline."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=-10.,
                     Sc=25.2,
                     rho_water=1025.)


def test_zero_wetted_area_hull_viscous():
    r"""Test a zero wetted area."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=0.,
                     rho_water=1025.)


def test_too_big_wetted_area_hull_viscous():
    r"""Test a too big wetted area."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=51.,
                     rho_water=1025.)


def test_wrong_multiplier_hull_viscous():
    r"""Test a wrong lwl multiplier."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     rho_water=1025.,
                     lwl_multiplier=0.69)
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     rho_water=1025.,
                     lwl_multiplier=1.01)


def test_negative_transition_rn_hull_viscous():
    r"""Negative transition reynolds."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     rho_water=1025.,
                     transition_reynolds_number=-100)


def test_wrong_formfactor_hull_viscous():
    r"""Test a wrong form factor."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     rho_water=1025.,
                     form_factor=0.89)
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     rho_water=1025.,
                     form_factor=1.61)


def test_wrong_rhowater_hull_viscous():
    r"""Test a wrong rho water."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     rho_water=900.)
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     rho_water=1100.)


def test_zero_kinematic_viscosity_hull_viscous():
    r"""Test a zerp kinematic viscosity."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     kinematic_viscosity=0.)


def test_wrong_friction_line_hull_viscous():
    r"""Friction line unknown."""
    with pytest.raises(ValueError):
        hull_viscous(boatspeed=3.5,
                     coe=(0.5, 0., -0.02),
                     lwl=10.,
                     Sc=25.2,
                     turbulent_friction_line=test_zero_kinematic_viscosity_hull_viscous)

# Wrong inputs ballast


def test_negative_length_ballast_viscous():
    r"""Negative length."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=-0.37,
                        bulb_diameter=0.037,
                        bulb_wetted_area=0.02)


def test_zero_diameter_ballast_viscous():
    r"""Zero diameter."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=0.37,
                        bulb_diameter=0.,
                        bulb_wetted_area=0.02)


def test_diameter_too_big_ballast_viscous():
    r"""Diameter too big."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=0.37,
                        bulb_diameter=0.2,
                        bulb_wetted_area=0.02)


def test_zero_wetted_area_ballast_viscous():
    r"""Zero wetted area."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=0.37,
                        bulb_diameter=0.037,
                        bulb_wetted_area=0.)


def test_negative_transition_rn_ballast_viscous():
    r"""Negative transition reynolds."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=0.37,
                        bulb_diameter=0.037,
                        bulb_wetted_area=0.02,
                        transition_reynolds_number=-100.)


def test_wrong_rho_water_ballast_viscous():
    r"""Rho water outside acceptable range."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=0.37,
                        bulb_diameter=0.037,
                        bulb_wetted_area=0.02,
                        rho_water=900)
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=0.37,
                        bulb_diameter=0.037,
                        bulb_wetted_area=0.02,
                        rho_water=1100)


def test_zero_kinematic_viscosity_ballast_viscous():
    r"""Zero kinematic viscosity."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=0.37,
                        bulb_diameter=0.037,
                        bulb_wetted_area=0.02,
                        kinematic_viscosity=0.)


def test_wrong_friction_line_ballast_viscous():
    r"""Friction line not a good value."""
    with pytest.raises(ValueError):
        ballast_viscous(boatspeed=1.,
                        heel_angle=0.,
                        coe=(0.5, 0., -0.2),
                        bulb_length=0.37,
                        bulb_diameter=0.037,
                        bulb_wetted_area=0.02,
                        turbulent_friction_line=test_zero_kinematic_viscosity_ballast_viscous)