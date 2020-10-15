#!/usr/bin/env python
# coding: utf-8

r"""Tests for the savitski.py module."""

from ydeos_hydrodynamics.savitsky import savitsky


def test_savitsky_happy_path():
    r"""Values from p204 of principles of yacht design"""
    boatspeed = 20.6
    m = 9750.
    LCG = 3.55
    VCG = 1.145
    b = 2.74
    epsilon = 3.
    beta = 17.9
    f = 0.49
    zero_bow_down_moment_trim_angle, M, Mh, Mf, Ma, R, Rf, Ra, ff, fa, Cv, \
    Cl_beta, Cl_zero, lambda_, delta_lambda, Lcp, e = \
        savitsky(boatspeed, m, LCG, VCG, b, epsilon, beta, f)

    assert zero_bow_down_moment_trim_angle == 4.168706856548404
    assert M == 0.029184050675667095
    assert Mh == -2796.0105530505853
    assert Mf == 2796.039737101261
    assert Ma == -0.000
    assert R == 13325.249564787735
    assert Rf == 6461.25223016241
    assert Ra == 0.000
    assert ff == 0.9237510284481013
    assert fa == 0.000
    assert Cv == 3.9740342082443276
    assert Cl_beta == 0.05855929710250599
    assert Cl_zero == 0.08508584811026254
    assert lambda_ == 1.8124177483423276
    assert delta_lambda == 0.000
    assert Lcp == 3.5435319901894853
    assert e == 0.006468009810514541
