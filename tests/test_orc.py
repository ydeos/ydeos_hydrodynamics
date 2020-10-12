#!/usr/bin/env python
# coding: utf-8

r"""Tests for the orc.py module"""

import pytest
from ydeos_hydrodynamics.orc import hull_residuary_resistance_orc_2013


# hrrorc2013 = HullResiduaryResistanceORC2013(Vc=7.63,
#                                             lwl=10.02,
#                                             bwl=3.17,
#                                             Tc=0.57)


def test_known_value_hull_residuary_resistance_orc_2013():
    r"""p 75 of Principles of Yacht Design - Larsson / Eliasson
    611 N @ Fn= 0.35 (lwl = 10.02 -> 0.35 * sqrt(9.81 * lwl) -> 3.47 m/s
    assertAlmostEquals(hrrorc2013.force(boatspeed = 3.47).force[0], -611.48676783980852)"""
    assert hull_residuary_resistance_orc_2013(boatspeed=3.47,
                                              Vc=7.63,
                                              lwl=10.02,
                                              bwl=3.17,
                                              Tc=0.57).fx == -522.2262497676658


def test_negative_lwl_hull_residuary_resistance_orc_2013():
    r"""Negative lwl"""
    with pytest.raises(ValueError):
        hull_residuary_resistance_orc_2013(boatspeed=3.47,
                                           Vc=7.63,
                                           lwl=-10.02,
                                           bwl=3.17,
                                           Tc=0.57)


def test_resistance_sign_and_symmetry_hull_residuary_resistance_orc_2013():
    r"""Test signs and symmetry of forces"""
    force = hull_residuary_resistance_orc_2013(boatspeed=2.,
                                               Vc=7.63,
                                               lwl=10.02,
                                               bwl=3.17,
                                               Tc=0.57)
    force1 = hull_residuary_resistance_orc_2013(boatspeed=-2.,
                                                Vc=7.63,
                                                lwl=10.02,
                                                bwl=3.17,
                                                Tc=0.57)
    assert force.fx == -force1.fx
    assert force.fx < 0