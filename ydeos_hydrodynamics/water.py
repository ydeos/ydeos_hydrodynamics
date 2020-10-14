# coding: utf-8

r"""Water characteristics."""

import numpy as np
from scipy.interpolate import Akima1DInterpolator

# rho water realistic limits
RHO_WATER_MIN = 950.
RHO_WATER_MAX = 1050.


def water_fresh_density(temperature_c: float) -> float:
    r"""Fresh water density.

    cf. https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html

    """
    t = np.array([0.1,     1.,     4.,    10.,    15.,    20.,    25.,    30.,    35.,    40.])
    d = np.array([999.85, 999.90, 999.97, 999.70, 999.10, 998.21, 997.05, 995.65, 994.03, 992.22])
    return Akima1DInterpolator(t, d)(temperature_c)


def water_sea_density(temperature_c: float) -> float:
    r"""For a salinity [g/lwl] of 35 g/lwl.

    https://www.engineeringtoolbox.com/sea-water-properties-d_840.html

    """
    t = np.array([0.,   10.,   15.,   20.,   23.,   26.,    30.])
    d = np.array([1028., 1027., 1026., 1025., 1024., 1023.,  1022.])
    return Akima1DInterpolator(t, d)(temperature_c)


def water_fresh_dynamic_viscosity(temperature_c: float) -> float:
    r"""Dynamic viscosity of fresh water (µ [Ns/m2])."""
    t = np.array([0.,       10.,       20.,       25.,       30.,       40.])
    d = np.array([0.0017914, 0.0013060, 0.0010016, 0.0008900, 0.0007972, 0.0006527])
    return Akima1DInterpolator(t, d)(temperature_c)


def water_sea_dynamic_viscosity(temperature_c: float) -> float:
    r"""Dynamic viscosity of sea water (µ [Ns/m2]).

    https://www.engineeringtoolbox.com/sea-water-properties-d_840.html
    salinity = 35g/lwl

    """
    t = np.array([0.,      5.,     10.,     15.,     20.,     25.,     30.])
    d = np.array([0.00188, 0.00162, 0.00141, 0.00123, 0.00108, 0.00096, 0.00086])
    return Akima1DInterpolator(t, d)(temperature_c)


def water_fresh_kinematic_viscosity(temperature_c: float) -> float:
    r"""Dynamic viscosity of fresh water [m2/s]."""
    return water_fresh_dynamic_viscosity(temperature_c) / water_fresh_density(temperature_c)


def water_sea_kinematic_viscosity(temperature_c: float) -> float:
    r"""Dynamic viscosity of sea water [m2/s]."""
    return water_sea_dynamic_viscosity(temperature_c) / water_sea_density(temperature_c)
