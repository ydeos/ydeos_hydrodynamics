#!/usr/bin/env python
# coding: utf-8

r"""Hypothesis tests for the savitsky.py module"""

from hypothesis import given
import hypothesis.strategies as st

from ydeos_hydrodynamics.savitsky import _cl_beta_cl_zero, _lambda


@given(boatspeed=st.floats(min_value=0.01, exclude_min=True, max_value=100.),
       m=st.floats(min_value=0.01, exclude_min=True, max_value=500000000.),
       b=st.floats(min_value=0.01, exclude_min=True, max_value=500.),
       beta=st.floats(min_value=0., max_value=90.))
def test_cl_beta_cl_zero(boatspeed, m, b, beta):
    cl_beta, cl_zero = _cl_beta_cl_zero(boatspeed, m, b, beta)
    assert cl_beta > 0.
    assert cl_zero > 0.


@given(trim_angle=st.floats(min_value=0., max_value=45.),
       Cv=st.floats(min_value=0.001, exclude_min=True, max_value=100.),
       Cl_zero=st.floats(min_value=0.001, exclude_min=True))
def test_lambda(trim_angle, Cv, Cl_zero):
    lambda_ = _lambda(trim_angle, Cv, Cl_zero)
    assert lambda_ > 0.
