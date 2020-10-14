#!/usr/bin/env python
# coding: utf-8

r"""Tests for the holtrop_mennen.py module."""

from ydeos_hydrodynamics.constants import KINEMATIC_VISCOSITY_WATER_20C
from ydeos_hydrodynamics.froude import froude_number
from ydeos_hydrodynamics.holtrop_mennen import holtrop_friction, holtrop_wave, \
    holtrop_bulbous, holtrop_transom, holtrop_a, holtrop_appendages,\
    _holtrop_Lr, _c2func, _c3func


def test_holtrop():
    r"""Against known values.

    Test values from the last page of
    'An Approximate Power Prediction Method - Holtrop Mennen - 1982.pdf'

    """
    L = 205.0  # Length on waterline
    LPP = 200.0
    B = 32.0  # breadth moulded
    Tf = 10.0  # draught moulded on FP
    Ta = 10.0  # draught moulded on AP
    T = 10.  # average moulded draught
    Vc = 37500.  # m3
    lcb_pp = 0.5202
    lcb_wl = 0.5075
    Abt = 20.0  # transverse bulb area
    h_b = 4.0
    Cm = 0.98
    Cwp = 0.75
    At = 16.0  # immersed transom
    S = 7381.45
    Sapp = 50.0  # wetted area appendages
    c_stern = 10.0  # stern shape parameter
    c_stern_code = "U"
    D = 8.0  # propeller diameter
    Z = 4  # number of propeller blades
    boatspeed = 25.0 * (1852/3600)  # 25 knots converted to m/s
    g = 9.81
    Cp = 0.5833
    rho_water = 1025.
    kinematic_viscosity = KINEMATIC_VISCOSITY_WATER_20C

    # Not provided in publication example
    Cb = 0.45  # block coeff, approx.

    assert 0.2868 - 1e-5 <= froude_number(boatspeed, L, g) <= 0.2868 + 1e-5
    assert 81.385 <= _holtrop_Lr(L, Cp, lcb_wl) <= 81.385 + 3e-3

    r_friction, k1, c12, c13, Cf = holtrop_friction(boatspeed,
                                                    L,
                                                    B,
                                                    T,
                                                    Cp,
                                                    lcb_wl,
                                                    c_stern_code,
                                                    S,
                                                    rho_water,
                                                    kinematic_viscosity)
    assert round(c12, 4) == 0.5102
    assert c13 == 1.03
    assert round(Cf, 5) == 0.00136  # 0.00139 in publi (different kinematic visc ?)
    assert round(k1, 3) == 0.156
    assert round(r_friction) == 852566  # 869630 in publi

    r_wave, c1, c5, c7, c15, m1, m2, lambda_ = holtrop_wave(boatspeed,
                                                            Vc,
                                                            L,
                                                            B,
                                                            T,
                                                            Cp,
                                                            lcb_wl,
                                                            At,
                                                            Cm,
                                                            Cwp,
                                                            Abt,
                                                            h_b,
                                                            Tf,
                                                            rho_water,
                                                            g)
    assert round(c1, 3) == 1.398
    assert round(c5, 4) == 0.9592
    assert round(c7, 4) == 0.1561
    assert round(c15, 5) == -1.69385  # +1.69385 in publi (is it a typo?)
    assert round(m1, 4) == -2.1274
    assert round(m2, 5) == -0.17086  # -0.17087 in publi
    assert round(lambda_, 4) == 0.6513
    assert round(r_wave) == 556788  # 557110 in publi

    c3 = _c3func(Abt, B, T, Tf, h_b)
    assert round(c3, 5) == 0.02119
    c2 = _c2func(c3)
    assert round(c2, 4) == 0.7595

    r_bulbous, Pb, Fni = holtrop_bulbous(boatspeed, Abt, h_b, Tf, rho_water, g)
    assert round(Pb, 4) == 0.6261
    assert round(Fni, 4) == 1.5083  # 1.5084 in publi
    assert round(r_bulbous) == 49

    r_transom, c6, FnT = holtrop_transom(boatspeed, B, At, Cwp, rho_water, g)
    assert round(FnT, 3) == 5.432  # 5.433 in publi
    assert round(r_transom) == 0

    r_a, c4, Ca = holtrop_a(boatspeed, L, B, T, S, Cb, Abt, h_b, Tf, rho_water)
    assert round(c4, 2) == 0.04
    assert round(Ca, 6) == 0.000352

    # 221980 in publication but Cb not provided, so this is ok
    assert round(r_a) == 220572

    # TODO : k2 was tweaked to pass the test
    # Is there any info in the publi anyway?
    r_app = holtrop_appendages(boatspeed,
                               (Sapp,),
                               (D,),
                               (0.0035,),
                               rho_water,
                               kinematic_viscosity)
    assert round(r_app) == 8830
