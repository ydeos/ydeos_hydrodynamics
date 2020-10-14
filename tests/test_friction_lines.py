# coding: utf-8

r"""Tests for the friction_lines.py module"""

import pytest

from ydeos_hydrodynamics.friction_lines import laminar, hughes, ittc57, \
    prandtl, appendage_friction_coefficient, bulb_friction_coefficient


def test_laminar():
    r"""Laminar friction line."""
    assert laminar(0., 10.) == 0.
    assert laminar(1., 10.) == 0.00042086558229325433


def test_ittc57():
    r"""ITTC57 turbulent friction line."""
    assert ittc57(0., 10.) == 0.
    assert ittc57(1., 10.) == 0.0030020815379449236


def test_hughes():
    r"""Hughes turbulent friction line."""
    assert hughes(0., 10.) == 0.
    assert hughes(1., 10.) == 0.002673832548653033


def test_prandtl():
    r"""Prandtl transitional friction line."""
    assert prandtl(1., 10., transition_reynolds_number=0.,
                   turbulent_friction_line=ittc57) == ittc57(1., 10.)
    assert prandtl(1., 10., transition_reynolds_number=0.,
                   turbulent_friction_line=hughes) == hughes(1., 10.)
    # transition never happens
    assert prandtl(1., 10., transition_reynolds_number=1e20,
                   turbulent_friction_line=ittc57) == laminar(1., 10.)

    with pytest.raises(ValueError):
        prandtl(1., 0., transition_reynolds_number=1e6,
                turbulent_friction_line=ittc57)
    with pytest.raises(ValueError):
        prandtl(1., 10., transition_reynolds_number=-1,
                turbulent_friction_line=ittc57)
    with pytest.raises(ValueError):
        prandtl(1., 10., transition_reynolds_number=1e6,
                turbulent_friction_line=test_hughes)


def test_appendage_cf():
    r"""ORC appendage Cf"""
    assert appendage_friction_coefficient(1., 10., thickness_to_chord=0.1) == 0.003531998254362007

    # foil too thick
    with pytest.raises(ValueError):
        appendage_friction_coefficient(1., 10., thickness_to_chord=0.3)

    # foil too thin
    with pytest.raises(ValueError):
        appendage_friction_coefficient(1., 10., thickness_to_chord=0.02)


def test_ballast_cf():
    r"""ORC appendage Cf"""
    assert bulb_friction_coefficient(1., 10.) == 0.004339311313727978
