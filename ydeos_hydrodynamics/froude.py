# coding: utf-8

r"""Froude number computations"""

from math import sqrt
from ydeos_hydrodynamics.constants import GRAVITY_STANDARD


def froude_number(speed: float, lwl: float, gravity: float = GRAVITY_STANDARD) -> float:
    """Froude number from speed [m/s] and waterline length [m]"""
    # lwl <= 0 or gravity <=0 will be caught by sqrt ValueError or ZeroDivision
    return abs(speed) / sqrt(lwl * gravity)


def speed_ms(froude: float, lwl: float, gravity: float = GRAVITY_STANDARD) -> float:
    """Speed [m/s] from Froude number and waterline length [m]"""
    # negative lwl or gravity will raise because of sqrt
    return abs(froude) * sqrt(lwl * gravity)
