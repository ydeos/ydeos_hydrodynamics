#!/usr/bin/env python
# coding: utf-8

r"""Holtrop example."""

from ydeos_hydrodynamics.froude import froude_number
from ydeos_hydrodynamics.holtrop_mennen import holtrop_mennen

v = 25  # knots
fn = froude_number(speed=v*(1852/3600), lwl=205.)
(r_friction, one_plus_k1, r_appendages, r_wave, r_bulbous, r_transom, r_a), total = \
    holtrop_mennen(v * (1852 / 3600),  # m/s
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
                   k2_app=(.5,))

print(f"        Froude number : {fn:.6f}")
print("\nResistance values in newtons [N]\n")
print(f"           R friction : {r_friction:.3f} [{r_friction / total * 100:.2f} %]")
print(f"               1 + k1 : {one_plus_k1}")
print(f"(1 + k1) * R friction : {one_plus_k1 * r_friction:.3f} [{one_plus_k1 * r_friction / total * 100:.2f} %]")
print(f"         R appendages : {r_appendages:.3f} [{r_appendages / total * 100:.2f} %]")
print(f"               R wave : {r_wave:.3f} [{r_wave / total * 100:.2f} %]")
print(f"            R bulbous : {r_bulbous:.3f} [{r_bulbous / total * 100:.2f} %]")
print(f"            R transom : {r_transom:.3f} [{r_transom / total * 100:.2f} %]")
print(f"      R a (roughness) : {r_a:.3f} [{r_a / total * 100:.2f} %]")
print(f"              R total : {total:.3f}")
