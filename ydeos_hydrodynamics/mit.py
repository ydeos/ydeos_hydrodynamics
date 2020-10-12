# coding: utf-8

r"""MIT resistance estimates"""

from typing import Union, Tuple
from math import sqrt
import numpy as np
from scipy.interpolate import interp1d
from ydeos_hydrodynamics.memoize import memoize
from ydeos_hydrodynamics.constants import GRAVITY_STANDARD, RHO_SEA_WATER_20C
from ydeos_hydrodynamics.water import RHO_WATER_MIN, RHO_WATER_MAX
from ydeos_hydrodynamics.force import Force
from ydeos_hydrodynamics.froude import froude_number


@memoize
def _get_hull_residuary_resistance_mit_interpolator(Vc: float,
                                                    lwl: float,
                                                    bwl: float,
                                                    Tc: float,
                                                    rho_water: float) -> interp1d:

    froude_numbers_table = [0.00, 0.1370668, 0.2056002, 0.2741336, 0.342667, 0.4112004, 0.4797338, 0.5482672, 0.6168006]

    a1 = [0., 9.8, 29.1, 47.2, 117.9, 413.7, 2427.6, 7421.8, 14055.]
    a2 = [0., 0.414, 0.617, 0.947, 0.831, 0.506, 0.085, 0.046, 0.032]
    a3 = [0., 0., 0., 0., 0., 0., 9., 49., 69.]

    Rrh = np.zeros(9)

    # Cv = Vc / lwl**3 * 1e-1
    Cv = 34977 * Vc * lwl**(-3)

    for i, froude in enumerate(froude_numbers_table):
        Rrh[i] = Vc * rho_water * GRAVITY_STANDARD * 1e-5 * a1[i] * (bwl/Tc)**a2[i] * Cv / sqrt(Cv**2 + a3[i])

    return interp1d(froude_numbers_table, Rrh, kind='slinear', bounds_error=True)


def hull_residuary_resistance_mit(boatspeed: float,
                                  Vc: float,
                                  lwl: float,
                                  bwl: float,
                                  Tc: float,
                                  coe: Union[Tuple[float, float, float], None] = None,
                                  rho_water: float = RHO_SEA_WATER_20C) -> Force:
    """The hull upright residuary resistance [N] according to the MIT regression

    boatspeed : boat speed [m/s]
    Vc : Canoe body volume of displaced water [m**3] , must be > 0.
    lwl : Length at Waterline [m] , must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    Tc : Canoe body draft [m], must be >= 0
    coe : tuple of 3 floats, optional (default [lwl / 2, 0, -Tc / 3]).
        Represents the x, y, z coordinates of the centre of effort.
        For residuary resistance, the surfacic barycentre of the underwater hull
        is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    rho_water : Water density [kg/m**3], must be >= 0.

    Returns a Force object, representing the residuary resistance [N] on the hull.

    """
    if Vc <= 0:
        raise ValueError("Vc should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    coe = (lwl / 2., 0., -Tc / 3.) if coe is None else coe

    interp = _get_hull_residuary_resistance_mit_interpolator(Vc, lwl, bwl, Tc, rho_water)

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.
    froude = froude_number(speed=boatspeed, lwl=lwl)
    if froude <= 0.6168006:
        residuary_resistance = max([0., interp(froude)])
    else:
        at_0_6068006 = interp(0.6068006)
        at_0_6168006 = interp(0.6168006)
        residuary_resistance = (at_0_6168006 +
                                (froude - 0.6168006) *
                                (at_0_6168006 - at_0_6068006) / 0.01)

    return Force(-residuary_resistance * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])
