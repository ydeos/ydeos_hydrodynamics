#!/usr/bin/env python
# coding: utf-8

r"""friction examples"""

from ydeos_hydrodynamics.friction_lines import prandtl, ittc57, hughes, \
    laminar, appendage_friction_coefficient, bulb_friction_coefficient

l = 2
s = 1
print("Cf (prandtl) fully turbulent : %f" % prandtl(speed=s,
                                                    dimension=l,
                                                    transition_reynolds_number=0))
print("                   Cf ittc57 : %f" % ittc57(speed=s, dimension=l))
print("  Cf (prandtl) fully laminar : %f" % prandtl(speed=s,
                                                    dimension=l,
                                                    transition_reynolds_number=1e12))
print("                  Cf laminar : %f" % laminar(speed=s, dimension=l))
print("                   Cf hughes : %f" % hughes(speed=s, dimension=l))
print("                Cf appendage : %f" %
      appendage_friction_coefficient(speed=s,
                                     dimension=l,
                                     thickness_to_chord=0.1))

print("                  Cf ballast : %f" % bulb_friction_coefficient(speed=s,
                                                                      dimension=l))
