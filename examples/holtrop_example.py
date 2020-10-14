#!/usr/bin/env python
# coding: utf-8

r"""Holtrop example."""

from ydeos_hydrodynamics.froude import froude_number
from ydeos_hydrodynamics.holtrop_mennen import holtrop_mennen

v = 25  # knots
print(froude_number(speed=v*(1852/3600), lwl=205.))
print(holtrop_mennen(v*(1852/3600),  # m/s
                     Vc=37500,
                     lwl=205.,
                     B=32.,
                     T=10.,
                     Cp=0.5833,
                     lcb_fraction=0.5075,
                     afterbody_shape="U",
                     S=7381.45,
                     At=16.,
                     Cm=0.98,
                     Cb=0.45,
                     Cwp=0.75,
                     Abt=20.,
                     h_B=4.0,
                     Tf=10.,
                     S_app=(50.,),
                     d_app=(200.,),
                     k2_app=(.5,)))
