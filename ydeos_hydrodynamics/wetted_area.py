# coding: utf-8

r"""Wetted area estimates"""

from scipy.interpolate import interp1d
from ydeos_hydrodynamics.memoize import memoize


@memoize
def _get_interpolator_wetted_area(upright_wsa: float, bwl: float, Tc: float, Cm: float) -> interp1d:
    heel_angles = [0.00, 5.00, 10.00, 15.00, 20.00, 25.00, 30.00, 35.00]

    s0 = [0.000, -4.112, -4.522, -3.291, 1.850, 6.510, 12.334, 14.648]
    s1 = [0.000, 0.054, -0.132, -0.389, -1.200, -2.305, -3.911, -5.182]
    s2 = [0.000, -0.027, -0.077, -0.118, -0.109, -0.066, 0.024, 0.102]
    s3 = [0.000, 6.329, 8.738, 8.949, 5.364, 3.443, 1.767, 3.497]

    wsa_heeled = [upright_wsa *
                  (1 +
                   0.01 * (s0[i] +
                           s1[i] * bwl / Tc +
                           s2[i] * (bwl / Tc) ** 2 +
                           s3[i] * Cm))
                  for i in range(8)]

    return interp1d(heel_angles, wsa_heeled, kind='cubic', bounds_error=False,
                    fill_value=wsa_heeled[-1])


def wetted_area(upright_wsa: float, bwl: float, Tc: float, Cm: float, heel: float) -> float:
    r"""Estimated wetted area for a given heel_angle
    Estimation of wetted surface area when heeled
    Value out of -35° <-> 35° of heel is equal to the value at 35°

    upright_wsa : Wetted area when the boat is upright
    bwl : Beam at waterline
    Tc : Canoe body draft
    Cm : Midship coefficient
    heel : Heel angle in degrees

    """
    if heel < 0:
        heel = -heel
    interp = _get_interpolator_wetted_area(upright_wsa, bwl, Tc, Cm)
    return interp(heel)


# the class version should now be slower because of the overhead than the function that calls a
# memoized _get_interpolator -> delete class
#
#
# class WettedArea(object):
#     """Faster implementation of the wetted_area() function.
#     The interpolator is only created once for a set of hull
#     characteristics
#
#     upright_wsa : Upright wetted area
#     bwl : Beam at waterline
#     Tc : Canoe body draft
#     Cm : Midship coefficient
#
#     """
#     def __init__(self, upright_wsa: float, bwl: float, Tc: float, Cm: float):
#         self.interp = _get_interpolator_wetted_area(upright_wsa, bwl, Tc, Cm)
#
#     def wetted_area(self, heel: float) -> float:
#         r"""Estimated wetted area at specified heel angle
#         heel : The heel angle in degrees
#
#         """
#         if heel < 0.:
#             heel = -heel
#         return self.interp(heel)
