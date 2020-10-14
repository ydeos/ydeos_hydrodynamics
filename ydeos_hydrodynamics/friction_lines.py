# coding: utf-8

r"""Friction lines."""

from typing import Callable
from math import sqrt, log10
import numpy as np
from scipy.interpolate import RectBivariateSpline, UnivariateSpline
from ydeos_hydrodynamics.reynolds import reynolds_number, transition_x
from ydeos_hydrodynamics.constants import KINEMATIC_VISCOSITY_WATER_20C


def laminar(speed: float,
            dimension: float,
            kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> float:
    """Cf (Friction Coefficient) for a laminar flow.

    The formula is known as the Blasius formula,
    and suffers very little controversy.
    speed : The flow speed [m/s]
    dimension : Characteristic dimension [m]
    kinematic_viscosity defaults to  water at 20°C : 1,004*10E-6 [m**2/s]

    """
    rn = reynolds_number(speed, dimension, kinematic_viscosity)

    if rn == 0.:
        return 0.
    # 1.32824 is for global Cf,
    # 0.664 is for a local Cf that requires integration over the surface
    return 1.32824 / sqrt(rn)


def ittc57(speed: float,
           dimension: float,
           kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> float:
    """Friction Coefficient according to ITTC-57.

    A 1.00 form factor should be used with this friction line result.
    speed : The flow speed [m/s]
    dimension : Characteristic dimension [m]
    kinematic_viscosity defaults to water at 20°C : 1,004*10E-6 [m**2/s]

    """
    rn = reynolds_number(speed, dimension, kinematic_viscosity)
    cf = 0. if rn == 0. else 0.075 / (log10(rn) - 2) ** 2
    return cf


def hughes(speed: float,
           dimension: float,
           kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> float:
    """Friction Coefficient according to the Hughes friction line.

    The Hughes friction line replaced the ittc57 friction lin
     in the 2013 ORC VPP
    A 1.05 form factor should be used with this friction line result.
    speed : The flow speed [m/s]
    dimension : Characteristic dimension [m]
    kinematic_viscosity defaults to water at 20°C : 1,004*10E-6 [m**2/s]

    """
    rn = reynolds_number(speed, dimension, kinematic_viscosity)
    cf = 0. if rn == 0. else 0.066 / (log10(rn) - 2.03) ** 2
    return cf


def prandtl(speed: float,
            dimension: float,
            transition_reynolds_number: float = 5e5,
            kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
            turbulent_friction_line: Callable = ittc57) -> float:
    """Friction coefficient, Prandtl formula (transitional friction line).

    This is best adapted for a situation with laminar and turbulent
    flow on the same body.

    speed : The flow speed [m/s]
    dimension : Characteristic dimension [m]
    transition_reynolds_number : the Reynolds number at which the transition
                                 occurs (default 5e5)
    kinematic_viscosity defaults to water at 20°C : 1,004*10E-6 [m**2/s]
    turbulent_friction_line : the function used to compute the turbulent Cf
                              (ittc57, hughes)

    """
    if dimension <= 0:
        raise ValueError("dimension should be strictly positive")
    if transition_reynolds_number < 0:
        raise ValueError("transition_reynolds_number should be positive or zero")
    if turbulent_friction_line not in [ittc57, hughes]:
        raise ValueError("The turbulent friction line should be ittc57 or hughes")

    # find transition position
    x_transition = transition_x(speed,
                                transition_reynolds_number,
                                kinematic_viscosity)

    if x_transition >= dimension:
        return laminar(speed, dimension)

    # Calculate Cf_lam_Xtr, Cf_turb_Xtr, Cf_turb_l
    cf_lam_xtr = laminar(speed, x_transition, kinematic_viscosity)
    cf_turb_xtr = turbulent_friction_line(speed,
                                          x_transition,
                                          kinematic_viscosity)
    cf_turb_l = turbulent_friction_line(speed, dimension, kinematic_viscosity)
    # Prandtl's formula for laminar + turbulent flow
    return cf_turb_l + (cf_lam_xtr - cf_turb_xtr) * x_transition / dimension


def appendage_friction_coefficient(speed: float,
                                   dimension: float,
                                   thickness_to_chord: float,
                                   kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> float:
    """Friction coefficient for an appendage.

    See page 53 of 2013 ORC VPP.

    speed : The flow speed [m/s]
    dimension : Characteristic dimension [m]
    thickness_to_chord : Thickness to chord ratio of the appendage
    kinematic_viscosity : defaults to water at 20°C : 1,004*10E-6 [m**2/s]

    """
    if not 0.03 < thickness_to_chord < 0.2:
        raise ValueError("thickness_to_chord is unrealistic")

    rn = reynolds_number(speed, dimension, kinematic_viscosity)

    # TODO : the interpolator is created at each function call -> memoize
    reynolds_numbers = [3.162e3, 1.00e4, 3.162e4, 1.0e5, 3.162e5, 1.00e6, 2.512e6, 6.310e6, 1.585e7, 5.012e7, 1.995e8]
    thickness_to_chord_ratios = [0., 0.1, 0.2]

    flat_plate_cfs = [24.85e-3, 13.86e-3, 7.73e-3, 4.95e-3, 3.46e-3, 3.00e-3, 3.0e-3, 3.0e-3, 2.81e-3, 2.39e-3, 1.96e-3]
    tc01_cfs = [42.07e-3, 28.93e-3, 20.20e-3, 10.74e-3, 4.99e-3, 3.62e-3, 3.62e-3, 3.62e-3, 3.39e-3, 2.88e-3, 2.36e-3]
    tc02_cfs = [44.12e-3, 30.51e-3, 21.42e-3, 11.50e-3, 5.40e-3, 3.94e-3, 3.94e-3, 3.94e-3, 3.69e-3, 3.14e-3, 2.59e-3]

    data = [flat_plate_cfs, tc01_cfs, tc02_cfs]

    sp = RectBivariateSpline(np.array(reynolds_numbers),
                             np.array(thickness_to_chord_ratios),
                             np.array(data).T,
                             kx=1,
                             ky=1,
                             s=0)

    # x <-> Rn, y <-> t/c , z <-> Cf

    return sp(rn, thickness_to_chord)[0][0]


def bulb_friction_coefficient(speed: float,
                              dimension: float,
                              kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> float:
    """Friction coefficient for a bulb.

    See page 53 of 2013 ORC VPP

    speed : The flow speed [m/s]
    dimension : Characteristic dimension [m]
    kinematic_viscosity defaults to water at 20°C : 1,004*10E-6 [m**2/s]

    """
    rn = reynolds_number(speed, dimension, kinematic_viscosity)

    # TODO : the interpolator is created at each function call -> memoize
    reynolds_numbers = [3.162e3, 1.00e4, 3.162e4, 1.0e5, 3.162e5, 1.00e6, 2.512e6, 6.310e6, 1.585e7, 5.012e7, 1.995e8]
    bulb_cfs = [59.29e-3, 44.00e-3, 32.66e-3, 16.54e-3, 6.51e-3, 4.49e-3, 4.49e-3, 4.49e-3, 4.21e-3, 3.57e-3, 2.93e-3]

    # k and s are determined by trial and error to look like the linear interpolation
    interp = UnivariateSpline(reynolds_numbers, bulb_cfs, k=2, s=0.000001)
    # interp = interp1d(reynolds_numbers, bulb_Cfs, kind='slinear', bounds_error=False, fill_value=0.)

    return float(interp(rn))
