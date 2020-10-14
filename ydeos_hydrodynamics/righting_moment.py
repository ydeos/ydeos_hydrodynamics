# coding: utf-8

r"""Righting moment estimates."""

from typing import Tuple
from math import sin, radians
from ydeos_hydrodynamics.constants import RHO_SEA_WATER_20C, GRAVITY_STANDARD
from ydeos_hydrodynamics.force import Force
from ydeos_hydrodynamics.froude import froude_number
from ydeos_hydrodynamics.water import RHO_WATER_MIN, RHO_WATER_MAX


def rm_estimate_gerritsma(boatspeed: float,
                          heel_angle: float,
                          displacement: float,
                          lwl: float,
                          bwl: float,
                          Tc: float,
                          Gdwl: float,
                          rho_water: float = RHO_SEA_WATER_20C) -> Tuple[Force, Force]:
    """Heeling righting moment estimate, Delft.

    Estimate the righting moment [N.m] from hull shape parameters
    according to the logic published in 1993 by Gerritsma

    boatspeed : boat speed [m/s]
    heel_angle : The heel_angle [degrees], between -90 and 90.
        A negative value represents heel to windward
    displacement : The total yacht displacement [kg], must be > 0.
    lwl : The length at waterline [m], must be > 0.
    bwl : The beam at waterline [m], must be > 0.
    Tc : The canoe body draft [m], must be > 0
    Gdwl : The vertical centre of gravity position relative to datum waterline
        (negative if below waterline)
    rho_water : Water density [kg/m**3],
                must be between RHO_WATER_MIN and RHO_WATER_MAX.

    Returns a list of 2 Force objects, of equal magnitude but of vertical
    and opposite directions representing the righting moment.

    References
    ----------
    http://www.boatdesign.net/forums/attachments/sailboats/
    51748d1293586862-dr-gerritsma-righting-moment-dr.-
                                                  gerritsmas-righting-moment.jpg

    http://www.boatdesign.net/forums/sailboats/
                                         dr-gerritsma-righting-moment-36057.html

    """
    if displacement <= 0:
        raise ValueError("displacement should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between "
                         f"{RHO_WATER_MIN} and {RHO_WATER_MAX}")

    D2 = -0.0406 + 0.0109 * (bwl / Tc) - 0.00105 * (bwl / Tc) ** 2
    D3 = 0.0636 - 0.0196 * (bwl / Tc)

    heel_sign = heel_angle / abs(heel_angle) if heel_angle != 0. else 0.

    froude = froude_number(speed=boatspeed, lwl=lwl)

    # Righting moment at 0 speed (rest)
    static_rm_estimate = (displacement *
                          GRAVITY_STANDARD *
                          ((0.664 * Tc + 0.111 * bwl ** 2 / Tc - Gdwl)
                           * sin(radians(abs(heel_angle))) +
                           lwl * D3 * (radians(abs(heel_angle))) ** 2))

    # Righting moment modification due to forward motion
    dynamic_rm_estimate = (displacement *
                           GRAVITY_STANDARD *
                           (lwl *
                            (D2 * (radians(abs(heel_angle))) * froude)))

    rm_estimate = static_rm_estimate + dynamic_rm_estimate
    # return ForceAndMoment([0.,0.,0.],
    #                       [-rm_estimate * heel_sign,0.,0.],[0.,0.,0.])

    # Approximation of the point of application of the righting force:
    # X: 1/2 of lwl
    # Y: - UNIT LENGTH (gives the right moment when multiplied
    #                   by the Z component of the force)
    # Z: 0

    righting_force = Force(0., 0., -rm_estimate * heel_sign,
                           lwl / 2., -1., 0.)
    compensating_force = Force(0., 0., rm_estimate * heel_sign,
                               lwl / 2., 0, 0.,)

    return righting_force, compensating_force


def rm_orc2013(heel_angle: float,
               displacement: float,
               lwl: float,
               bwl: float,
               Tc: float,
               rho_water: float = RHO_SEA_WATER_20C) -> Tuple[Force, Force]:
    """Heeling righting moment estimate, ORC.

    Righting moment estimate, inspired by the formulation at page 30 of
    the 'ORC VPP Documentation 2013, updated july 4th'

    DSPM    Displacement in measurement trim [kg]
    IMSL    Effective sailing length [m] weighted average of lengths
            for 3 conditions of flotation
    BTR     Beam Depth Ratio [] - effective beam divided by
            the effective hull depth
    VOL     Canoe body volume [m**3]
    SA      Sail area upwind [m**2]
    HA      heeling arm [N*m], defined as :
            (CEH main*AREA main + CEH jib*AREA jib) / SA + HBI + DHKA*0.45
    HBI     Height [m] of base of I
    DHKA    Draft [m] of keel and hull adjusted
    B       Effective beam [m]
            @see www.orc.org/rules/ORC%20Rating%20Rules%202010.pdf (page 6)

    """
    if displacement <= 0:
        raise ValueError("displacement should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    a0 = -0.00410481
    a1 = -0.00003999
    a2 = -0.00017008
    a3 = 0.000019183
    a4 = 0.003602739

    heel_sign = heel_angle / abs(heel_angle) if heel_angle != 0. else 0.

    DSPM = displacement
    IMSL = lwl
    BTR = bwl / Tc
    VOL = displacement / rho_water
    B = bwl
    SA = 0
    HA = 0.7
    rm_default = 1.025 * DSPM * IMSL * (a0 +
                                        a1*BTR +
                                        a2*(VOL**(1./3.)/IMSL) +
                                        a3*SA*HA/B**3 + a4*B/VOL**(1./3.))

    righting_force = Force(0., 0., -rm_default * heel_sign, lwl / 2., -1., 0.)
    compensating_force = Force(0., 0., rm_default * heel_sign, lwl / 2., 0, 0.)

    return righting_force, compensating_force


# Trim righting moment estimate

# This is created by swapping beam and length in Gerritsma's heeling righting
# moment estimate


def trm_estimate(trim_angle: float,
                 displacement: float,
                 lwl: float,
                 bwl: float,
                 Tc: float,
                 Gdwl: float,
                 rho_water: float = RHO_SEA_WATER_20C) -> Tuple[Force, Force]:
    """Trim righting moment estimate.

    Estimate the righting moment [N.m] from hull shape parameters
    according to the logic published in 1993 by Gerritsma

    trim_angle : trim angle [degrees].
        Negative is bow down, Positive is bow up
    displacement : The total yacht displacement [kg], must be > 0.
    lwl : The length at waterline [m], must be > 0.
    bwl : The beam at waterline [m], must be > 0.
    Tc : The canoe body draft [m], must be > 0
    Gdwl : The vertical centre of gravity position relative to datum waterline
        (negative if below waterline)
    rho_water : Water density [kg/m**3],
                must be between RHO_WATER_MIN and RHO_WATER_MAX.

    Returns a list of 2 Force objects, of equal magnitude but of vertical
    and opposite directions representing the righting moment.

    References
    ----------
    http://www.boatdesign.net/forums/attachments/sailboats/
        51748d1293586862-dr-gerritsma-righting-moment-dr.-gerritsmas-righting-moment.jpg
    http://www.boatdesign.net/forums/sailboats/ dr-gerritsma-righting-moment-36057.html

    """
    if displacement <= 0:
        raise ValueError("displacement should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between "
                         f"{RHO_WATER_MIN} and {RHO_WATER_MAX}")

    D3 = 0.0636 - 0.0196 * (lwl / Tc)

    trim_sign = trim_angle / abs(trim_angle) if trim_angle != 0. else 0.

    # Righting moment at 0 speed (rest)
    static_trm_estimate = (displacement * GRAVITY_STANDARD *
                           ((0.664 * Tc + 0.111 * lwl ** 2 / Tc - Gdwl) *
                            sin(radians(abs(trim_angle))) + bwl * D3 * (radians(abs(trim_angle))) ** 2))

    trm_estimate_ = static_trm_estimate

    trim_righting_force = Force(0., 0., trm_estimate_ * trim_sign,
                                lwl / 2. - 1., 0., 0.)
    compensating_force = Force(0., 0., -trm_estimate_ * trim_sign,
                               lwl / 2., 0, 0.)

    return trim_righting_force, compensating_force
