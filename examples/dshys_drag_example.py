#!/usr/bin/env python
# coding: utf-8

r"""Sysser01 drag values from various regressions."""

import matplotlib.pyplot as plt
import numpy as np
from ydeos_hydrodynamics.delft import hull_residuary_resistance_ks_series4, \
    hull_residuary_resistance_la_series3
from ydeos_hydrodynamics.mit import hull_residuary_resistance_mit
from ydeos_hydrodynamics.orc import hull_residuary_resistance_orc_2013
from ydeos_hydrodynamics.froude import speed_ms, froude_number
from ydeos_hydrodynamics.viscous import hull_viscous

plot_type = "viscous+residuary"
# plot_type = "residuary"


class Sysser01:
    r"""Data container class for Sysser01."""
    Vc = 0.0376136
    lwl = 1.6
    bwl = 0.5072
    lcb_frac = 0.53664
    waterplane_area = 0.558566
    lcf_frac = 0.55328
    Tc = 0.12704
    wsa = 0.642534
    Cp = 0.5644135025


# fns = [0.2, 0.3, 0.4, 0.5]
fns = np.linspace(start=0., stop=0.8, num=100, endpoint=True)

sysser01_r_ks = list()
sysser01_r_la = list()
sysser01_r_orc = list()
sysser01_r_mit = list()

for fn in fns:
    print("Fn : %.2f" % fn)

    v = speed_ms(froude=fn, lwl=Sysser01.lwl)

    print("V : %.2f [m/s]" % v)

    sysser01_visc = -hull_viscous(v, Sysser01.lwl, Sysser01.wsa).fx

    sysser01_rrla = -hull_residuary_resistance_la_series3(
        boatspeed=v,
        Vc=Sysser01.Vc,
        lwl=Sysser01.lwl,
        bwl=Sysser01.bwl,
        Tc=Sysser01.Tc,
        Aw=Sysser01.waterplane_area,
        lcb_fraction=Sysser01.lcb_frac,
        Cp=Sysser01.Cp).fx

    sysser01_rrks = -hull_residuary_resistance_ks_series4(
        boatspeed=v,
        Vc=Sysser01.Vc,
        lwl=Sysser01.lwl,
        bwl=Sysser01.bwl,
        Aw=Sysser01.waterplane_area,
        Sc=Sysser01.wsa,
        lcb_fraction=Sysser01.lcb_frac,
        lcf_fraction=Sysser01.lcf_frac,
        Cp=Sysser01.Cp).fx

    sysser01_rrorc = -hull_residuary_resistance_orc_2013(boatspeed=v,
                                                         Vc=Sysser01.Vc,
                                                         lwl=Sysser01.lwl,
                                                         bwl=Sysser01.bwl,
                                                         Tc=Sysser01.Tc).fx

    sysser01_rrmit = -hull_residuary_resistance_mit(boatspeed=v,
                                                    Vc=Sysser01.Vc,
                                                    lwl=Sysser01.lwl,
                                                    bwl=Sysser01.bwl,
                                                    Tc=Sysser01.Tc).fx

    if plot_type == "viscous+residuary":
        sysser01_r_ks.append(sysser01_visc + sysser01_rrks)
        sysser01_r_la.append(sysser01_visc + sysser01_rrla)
        sysser01_r_orc.append(sysser01_visc + sysser01_rrorc)
        sysser01_r_mit.append(sysser01_visc + sysser01_rrmit)
    elif plot_type == "residuary":
        sysser01_r_ks.append(sysser01_rrks)
        sysser01_r_la.append(sysser01_rrla)
        sysser01_r_orc.append(sysser01_rrorc)
        sysser01_r_mit.append(sysser01_rrmit)
    else:
        raise ValueError("Unknown plot type")

    # print(60 * '*')
    # print("Fn : %f" % fn)
    # print("Speed : %f" % v)
    # print("")
    # print("Sysser01")
    # print(sysser01_visc)
    # print(sysser01_rrks)
    # print(sysser01_rrla)
    # print(sysser01_visc + sysser01_rrks)
    # print(sysser01_visc + sysser01_rrla)

plt.plot(fns, sysser01_r_ks,
         c="b",
         ls="--",
         label="Sysser01 viscous + residuary (ks)")
plt.plot(fns, sysser01_r_la,
         c="b",
         label="Sysser01 viscous + residuary (la)")

plt.plot(fns, sysser01_r_orc,
         c="orange",
         label="Sysser01 viscous + residuary (orc)")

plt.plot(fns, sysser01_r_mit,
         c="black",
         # ls="--",
         label="Sysser01 viscous + residuary (mit)")

if plot_type == "viscous+residuary":
    f = [froude_number(speed=v, lwl=Sysser01.lwl) for v in [0.396,
                                                            0.594,
                                                            0.792,
                                                            0.99,
                                                            1.189,
                                                            1.387,
                                                            1.585,
                                                            1.783,
                                                            1.981,
                                                            2.179,
                                                            2.377]]

    r = [0.257851, 0.616414, 1.1084, 1.88482, 3.01934, 4.75255, 9.79318,
         19.4122, 32.7027, 44.5576, 51.149]

    plt.scatter(f, r, label="Towing tank -Upright")

plt.legend()
plt.grid()
plt.title("R [N] - %s" % plot_type)
plt.xlabel("Fn")
plt.ylabel("Resistance")

plt.show()
