# coding: utf-8

r"""Geometric computations.

Todo:
----
coe moves around an arbitrary axis, moves with 3 angles ...
-> use transformations.py

"""

from typing import Tuple
from math import cos, sin, radians


def shift_coe_around_x(coe: Tuple[float, float, float],
                       heel_angle: float) -> Tuple[float, float, float]:
    r"""COE movement around x axis (heeling).

    Move the coe around the x axis to represent its displacement
    due to the heel angle

    """
    px, py, pz = coe
    return (px,
            py * cos(radians(heel_angle)) + pz * sin(radians(heel_angle)),
            pz * cos(radians(heel_angle)) - py * sin(radians(heel_angle)))
