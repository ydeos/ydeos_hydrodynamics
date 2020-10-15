#!/usr/bin/env python
# coding: utf-8

r"""Off limit behavior of the Delft regressions.

Create a 3d graph with x-> speed, y->bwl, z->wave resistance in order
to visualize the behaviour of the Delft series when out the prescribed limits

Todo:
----
improve the visualization

"""

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import numpy as np

from ydeos_hydrodynamics.delft import hull_residuary_resistance_ks_series4, \
    check_lwl_to_bwl_s4
from ydeos_hydrodynamics.mit import hull_residuary_resistance_mit
from ydeos_hydrodynamics.froude import speed_ms

# Base boat characteristics
lwl = 1.22
b = 0.18
t = 0.06
cm = 0.75
cp = 0.55

# Speeds in m/s, specified using the Froude numbers
speeds = np.linspace(speed_ms(0.1, lwl), speed_ms(0.5, lwl), 100)
# bwl varies from 0.1 to 0.3 meters
bwls = np.linspace(0.10, 0.3, 100)

base_volume = lwl * b * t * cm * cp

print(f"Base volume is {base_volume * 1000} [l]")

# data lists for data within limits
xdata, ydata, zdata = [], [],  []

# data lists for off limits data
xdata_off_limits, ydata_off_limits, zdata_off_limits = [], [], []

# mit wave resistance model for comparison
xdata_mit, ydata_mit, zdata_mit = [], [], []

for x in speeds:
    for bwl in bwls:
        r = -hull_residuary_resistance_ks_series4(boatspeed=x,
                                                  Vc=base_volume,
                                                  lwl=lwl,
                                                  bwl=bwl,
                                                  Aw=lwl * bwl * 0.7,
                                                  Sc=lwl * bwl * 0.8,
                                                  lcb_fraction=0.53,
                                                  lcf_fraction=0.56,
                                                  Cp=0.56,
                                                  coe=None,
                                                  rho_water=1025.).fx
        r_mit = -hull_residuary_resistance_mit(x,
                                               base_volume,
                                               lwl,
                                               bwl,
                                               t * b / bwl).fx
        xdata_mit.append(x)
        ydata_mit.append(bwl)
        zdata_mit.append(r_mit)

        if not check_lwl_to_bwl_s4(value=lwl / bwl):
            xdata_off_limits.append(x)
            ydata_off_limits.append(bwl)
            zdata_off_limits.append(r)
        else:
            xdata.append(x)
            ydata.append(bwl)
            zdata.append(r)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_trisurf(xdata,
                ydata,
                zdata,
                color='chartreuse',
                shade=False,
                label="within limits")
ax.plot_trisurf(xdata_off_limits,
                ydata_off_limits,
                zdata_off_limits,
                color='grey',
                shade=False,
                label="off limits")
ax.plot_trisurf(xdata_mit,
                ydata_mit,
                zdata_mit,
                color='darkblue',
                shade=False,
                label="mit")

ax.set_xlabel('Speed [m/s]')
ax.set_ylabel('Bwl [m]')
ax.set_zlabel('R [N]')

plt.show()
