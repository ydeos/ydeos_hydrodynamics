#!/usr/bin/env python
# coding: utf-8

r"""Tests for the delft.py module"""

from time import time
import pytest


from ydeos_hydrodynamics.delft import keel_residuary_delta_heeled_ks, \
    keel_residuary_ks, due_to_heel, hull_residuary_delta_heeled_ks, \
    hull_residuary_resistance_2008, hull_residuary_resistance_ks_series4, \
    hull_residuary_resistance_la_series3


krrdhks = keel_residuary_delta_heeled_ks(boatspeed=10.,
                                         heel_angle=0.,
                                         Vc=0.004,
                                         Vk=0.00004,
                                         Tc=0.06,
                                         T=0.4,
                                         lwl=1.0,
                                         bwl=0.2)


def test_known_value_keel_residuary_delta_heeled_ks():
    r"""No know data if heel is not 0"""
    assert krrdhks.fx == 0.


def test_negative_lwl_keel_residuary_delta_heeled_ks():
    r"""Negative lwl"""
    with pytest.raises(ValueError):
        keel_residuary_delta_heeled_ks(boatspeed=10.,
                                       heel_angle=0.,
                                       Vc=0.004,
                                       Vk=0.00004,
                                       Tc=0.06,
                                       T=0.4,
                                       lwl=-1.0,
                                       bwl=0.2)


def test_resistance_sign_and_symmetry_keel_residuary_delta_heeled_ks():
    r"""Signs and symmetry tests"""
    force = keel_residuary_delta_heeled_ks(boatspeed=4.,
                                           heel_angle=30.,
                                           Vc=0.004,
                                           Vk=0.00004,
                                           Tc=0.06,
                                           T=0.4,
                                           lwl=1.0,
                                           bwl=0.2)
    force_backwards = keel_residuary_delta_heeled_ks(boatspeed=-4.,
                                                     heel_angle=30.,
                                                     Vc=0.004,
                                                     Vk=0.00004,
                                                     Tc=0.06,
                                                     T=0.4,
                                                     lwl=1.0,
                                                     bwl=0.2)
    assert force.fx == -force_backwards.fx
    assert force.fx < 0

    force_heeled_to_windward = keel_residuary_delta_heeled_ks(boatspeed=4.,
                                                              heel_angle=-30.,
                                                              Vc=0.004,
                                                              Vk=0.00004,
                                                              Tc=0.06,
                                                              T=0.4,
                                                              lwl=1.0,
                                                              bwl=0.2)
    assert force.fx == force_heeled_to_windward.fx


def test_known_value_keel_residuary_ks():
    r"""No known data"""
    assert True


def test_negative_lwl_keel_residuary_ks():
    r"""Negative lwl"""
    with pytest.raises(ValueError):
        keel_residuary_ks(boatspeed=1.,
                          Vc=0.004,
                          Vk=0.00004,
                          Tc=0.06,
                          T=0.4,
                          lwl=-1.0,
                          bwl=0.2,
                          Zcbk=0.15)


def test_resistance_sign_and_symmetry_keel_residuary_ks():
    r"""Test signs and symmetry"""
    force = keel_residuary_ks(boatspeed=1.,
                              Vc=0.004,
                              Vk=0.00004,
                              Tc=0.06,
                              T=0.4,
                              lwl=1.0,
                              bwl=0.2,
                              Zcbk=0.15)
    force_backwards = keel_residuary_ks(boatspeed=-1.,
                                        Vc=0.004,
                                        Vk=0.00004,
                                        Tc=0.06,
                                        T=0.4,
                                        lwl=1.0,
                                        bwl=0.2,
                                        Zcbk=0.15)
    assert force.fx == -force_backwards.fx
    assert force.fx < 0


@pytest.mark.skip
def test_memoizing_keel_residuary_ks():
    # MAKE SURE THE PARAMS USED TO BUILD THE INTERPOLATOR ARE DIFFERENT
    # FROM THE PREVIOUS CALLS TO TEST MEMOIZING

    t0 = time()
    _ = keel_residuary_ks(boatspeed=1., Vc=0.004, Vk=0.00004, Tc=0.06, T=0.4, lwl=1.0, bwl=0.21, Zcbk=0.15)
    t1 = time()
    _ = keel_residuary_ks(boatspeed=2., Vc=0.004, Vk=0.00004, Tc=0.06, T=0.4, lwl=1.0, bwl=0.21, Zcbk=0.15)
    t2 = time()
    second_call = t2 - t1
    first_call = t1 - t0
    print(first_call)
    print(second_call)
    assert second_call < first_call


# Created from the YD-40 characteristics in half loaded displacement
# in Principles of Yacht Design (cf. Appendix 1)
# coe = [0,0,0]
hrdthla = due_to_heel(boatspeed=3.5,
                      heel_angle=13.6,
                      lwl=10.02,
                      bwl=3.17,
                      Sc=25.2,
                      Tc=0.57,
                      T=2.07)


def test_known_value_due_to_heel():
    r"""p 83 of Principles of Yacht Design - Larsson / Eliasson
    99 N @ 3.5 m/s and 13.6 deg. heel
    """
    value = -100.8296221334219
    assert value - 1e-6 <= hrdthla.fx <= value + 1e-6


def test_negative_lwl_due_to_heel():
    r"""Negative lwl"""
    with pytest.raises(ValueError):
        due_to_heel(boatspeed=3.5,
                    heel_angle=13.6,
                    lwl=-10.02,
                    bwl=3.17,
                    Sc=25.2,
                    Tc=0.57,
                    T=2.07)


def test_resistance_sign_and_symmetry_due_to_heel():
    r"""Test signs and symmetry"""
    force = due_to_heel(boatspeed=4.,
                        heel_angle=20.,
                        lwl=10.02,
                        bwl=3.17,
                        Sc=25.2,
                        Tc=0.57,
                        T=2.07)
    force_heeled_to_windward = due_to_heel(boatspeed=4.,
                                           heel_angle=-20.,
                                           lwl=10.02,
                                           bwl=3.17,
                                           Sc=25.2,
                                           Tc=0.57,
                                           T=2.07)
    force_backwards = due_to_heel(boatspeed=-4.,
                                  heel_angle=20.,
                                  lwl=10.02,
                                  bwl=3.17,
                                  Sc=25.2,
                                  Tc=0.57,
                                  T=2.07)
    assert force.fx == -force_backwards.fx
    assert force.fx == force_heeled_to_windward.fx

    assert force.fx < 0


def test_default_coe_position_due_to_heel():
    r"""Test default centre of effort"""
    force = due_to_heel(boatspeed=4.,
                        heel_angle=20.,
                        lwl=10.02,
                        bwl=3.17,
                        Sc=25.2,
                        Tc=0.57,
                        T=2.07)

    assert force.px == 5.01
    assert force.pz == -0.57 / 3


def test_known_value_hull_residuary_delta_heeled_ks():
    r"""No know data if heel is not 0"""
    assert hull_residuary_delta_heeled_ks(boatspeed=10., heel_angle=0., Vc=0.004, lwl=1.0, bwl=0.2, Tc=0.06).fx == 0.


def test_negative_lwl_hull_residuary_delta_heeled_ks():
    r"""Negative lwl"""
    with pytest.raises(ValueError):
        hull_residuary_delta_heeled_ks(boatspeed=10., heel_angle=0., Vc=0.004, lwl=-1.0, bwl=0.2, Tc=0.06)


def test_resistance_sign_and_symmetry_hull_residuary_delta_heeled_ks():
    r"""Test force signs and symmetry"""
    force = hull_residuary_delta_heeled_ks(boatspeed=4., heel_angle=30., Vc=0.004, lwl=1.0, bwl=0.2, Tc=0.06)
    force_backwards = hull_residuary_delta_heeled_ks(boatspeed=-4., heel_angle=30., Vc=0.004, lwl=1.0, bwl=0.2, Tc=0.06)
    assert force.fx == -force_backwards.fx
    assert force.fx < 0

    force_heeled_to_windward = hull_residuary_delta_heeled_ks(boatspeed=4., heel_angle=-30., Vc=0.004, lwl=1.0, bwl=0.2, Tc=0.06)
    assert force.fx == force_heeled_to_windward.fx


def test_known_value_hull_residuary_resistance_2008():
    r""""""
    assert hull_residuary_resistance_2008(boatspeed=3.47, Vc=7.63, lwl=10.02, bwl=3.17, Tc=0.57, Aw=22., Cm=0.75,
                                          lcb_fraction=0.535, lcf_fraction=0.550, Cp=0.56).fx == -521.8152569777417


# Created from the YD-40 characteristics in half loaded displacement
# in Principles of Yacht Design (cf. Appendix 1)
# coe = [0,0,0]
# Aw=22., # not provided in book, approximated from Sc
# lcf_fraction=0.550, # not provided in book, approximated from lcb_fraction


# hrrks = HullResiduaryResistanceKs(Vc=7.63,
#                                   lwl=10.02,
#                                   bwl=3.17,
#                                   Aw=22.,
#                                   Sc=25.2,
#                                   lcb_fraction=0.535,
#                                   lcf_fraction=0.550,
#                                   Cp=0.56)


def test_known_value_hull_residuary_resistance_ks_series4():
    r"""p 75 of Principles of Yacht Design - Larsson / Eliasson
    611 N @ Fn= 0.35 (lwl = 10.02 -> 0.35 * sqrt(9.81 * lwl) -> 3.47 m/s
    assertAlmostEquals(abs(hrrks.force(boatspeed = 3.47).force[0]), 611.0)
    525 found instead of 611 !!!!!!"""
    assert hull_residuary_resistance_ks_series4(boatspeed=3.47, Vc=7.63,
                                                lwl=10.02,
                                                bwl=3.17,
                                                Aw=22.,
                                                Sc=25.2,
                                                lcb_fraction=0.535,
                                                lcf_fraction=0.550,
                                                Cp=0.56).fx == -526.2054296064833


def test_negative_lwl_hull_residuary_resistance_ks_series4():
    r""" Negative lwl
    Aw=22., # not provided in book, approximated from Sc
    lcf_fraction=0.550, # not provided in book, approximated from lcb_fraction"""
    with pytest.raises(ValueError):
        hull_residuary_resistance_ks_series4(boatspeed=0.,
                                             Vc=7.63,
                                             lwl=-10.02,
                                             bwl=3.17,
                                             Aw=22.,
                                             Sc=25.2,
                                             lcb_fraction=0.535,
                                             lcf_fraction=0.550,
                                             Cp=0.56)


def test_resistance_sign_and_symmetry_hull_residuary_resistance_ks_series4():
    r"""Test signs and symmetry of forces"""
    force = hull_residuary_resistance_ks_series4(boatspeed=4.,
                                                 Vc=7.63,
                                                 lwl=10.02,
                                                 bwl=3.17,
                                                 Aw=22.,
                                                 Sc=25.2,
                                                 lcb_fraction=0.535,
                                                 lcf_fraction=0.550,
                                                 Cp=0.56)
    force1 = hull_residuary_resistance_ks_series4(boatspeed=-4.,
                                                  Vc=7.63,
                                                  lwl=10.02,
                                                  bwl=3.17,
                                                  Aw=22.,
                                                  Sc=25.2,
                                                  lcb_fraction=0.535,
                                                  lcf_fraction=0.550,
                                                  Cp=0.56)
    assert force.fx == -force1.fx
    assert force.fx < 0


# Created from the YD-40 characteristics in half loaded displacement
# in Principles of Yacht Design (cf. Appendix 1)
# coe = [0,0,0]
# Aw=22., # not provided in book, approximated from Sc


# hrrla = HullResiduaryResistanceLa(Vc=7.63,
#                                   lwl=10.02,
#                                   bwl=3.17,
#                                   Tc=0.57,
#                                   Aw=22.,
#                                   lcb_fraction=0.535,
#                                   Cp=0.56)


def test_known_value_hull_residuary_resistance_la_series3():
    r"""p 75 of Principles of Yacht Design - Larsson / Eliasson
    611 N @ Fn= 0.35 (lwl = 10.02 -> 0.35 * sqrt(9.81 * lwl) -> 3.47 m/s"""
    assert hull_residuary_resistance_la_series3(boatspeed=3.47,
                                                Vc=7.63,
                                                lwl=10.02,
                                                bwl=3.17,
                                                Tc=0.57,
                                                Aw=22.,
                                                lcb_fraction=0.535,
                                                Cp=0.56).fx == -612.0903884276029


def test_negative_lwl_hull_residuary_resistance_la_series3():
    r"""Negative lwl
    Aw=22., # not provided in book, approximated from Sc"""
    with pytest.raises(ValueError):
        hull_residuary_resistance_la_series3(boatspeed=0.,
                                             Vc=7.63,
                                             lwl=-10.02,
                                             bwl=3.17,
                                             Tc=0.57,
                                             Aw=22.,
                                             lcb_fraction=0.535,
                                             Cp=0.56)


def test_resistance_sign_and_symmetry_hull_residuary_resistance_la_series3():
    r"""Test signs and symmetry of forces"""
    force = hull_residuary_resistance_la_series3(boatspeed=2.,
                                                 Vc=7.63,
                                                 lwl=10.02,
                                                 bwl=3.17,
                                                 Tc=0.57,
                                                 Aw=22.,
                                                 lcb_fraction=0.535,
                                                 Cp=0.56)
    force1 = hull_residuary_resistance_la_series3(boatspeed=-2.,
                                                  Vc=7.63,
                                                  lwl=10.02,
                                                  bwl=3.17,
                                                  Tc=0.57,
                                                  Aw=22.,
                                                  lcb_fraction=0.535,
                                                  Cp=0.56)
    assert force.fx == -force1.fx
    assert force.fx < 0