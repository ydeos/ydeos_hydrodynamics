#!/usr/bin/env python
# coding: utf-8

r"""Friction lines examples."""

from ydeos_hydrodynamics.friction_lines import prandtl, ittc57, hughes, \
    laminar, appendage_friction_coefficient, bulb_friction_coefficient

length = 2
s = 1
print(f"Cf (prandtl) fully turbulent : "
      f"{prandtl(speed=s, dimension=length, transition_reynolds_number=0)}")
print(f"                   Cf ittc57 : {ittc57(speed=s, dimension=length)}")
print(f"  Cf (prandtl) fully laminar : "
      f"{prandtl(speed=s, dimension=length, transition_reynolds_number=1e12)}")
print(f"                  Cf laminar : {laminar(speed=s, dimension=length)}")
print(f"                   Cf hughes : {hughes(speed=s, dimension=length)}")
print(f"                Cf appendage : "
      f"{appendage_friction_coefficient(speed=s, dimension=length, thickness_to_chord=0.1)}")

print(f"                  Cf ballast : "
      f"{bulb_friction_coefficient(speed=s, dimension=length)}")
