# coding: utf-8

r"""Sideforce and corresponding resistance estimates."""

from typing import Tuple, Optional
from math import cos, radians, pi, sin, degrees, sqrt, tan

from ydeos_hydrodynamics.constants import RHO_SEA_WATER_20C
from ydeos_hydrodynamics.force import Force
from ydeos_hydrodynamics.water import RHO_WATER_MIN, RHO_WATER_MAX
from ydeos_hydrodynamics.geometry import shift_coe_around_x


def sideforce_by_element_fossati(boatspeed: float,
                                 heel_angle: float,
                                 leeway_angle: float,
                                 rudder_angle: float,
                                 coe: Tuple[float, float, float],
                                 projected_area: float,
                                 span: float,
                                 aspect_ratio_multiplier: float = 2.,
                                 use_cos_heel_angle: bool = True,
                                 rudder_angle_influence: float = 0.,
                                 speed_multiplier: float = 1.0,
                                 rho_water: float = RHO_SEA_WATER_20C,
                                 # for real 2D sections
                                 filename_2d_section: str = None,
                                 ncrit: float = 1.0,
                                 add_viscous_drag: bool = False) -> Force:
    """Lift (i.e. sideforce) and drag due to lift [N] of an underwater element,
    mostly due to leeway.
    Viscous effects (friction and pressure drag) are not accounted for.

    boatspeed : boat speed [m/s]
    heel_angle : The heel_angle [degrees], between -90 and 90.
    leeway_angle : The leeway angle
    rudder_angle : The rudder angle

    coe : tuple of 3 floats
        Represents the x, y, z coordinates of the centre of effort.
        The centre of effort can be estimated by projecting the projected
        surface barycentre to the 1/4 chord line.
    projected_area : Projected area [m**2] of the underwater element,
        must be >= 0
    span : Span [m] of the underwater element, must be >= 0.
    aspect_ratio_multiplier : Modifier of the geometric aspect ratio calculated
        from projected area and span to get effective aspect ratio
        Must be between 0.5 and 2.5
        Use 2. for a lifting plane with a plate effect at one end.
    use_cos_heel_angle : If True, multiply the angle of attack by the cosine of the heel angle.
        If False, use the raw angle of attack.
        Use False for the hull / canoe body; use True for keel and rudder.
    rudder_angle_influence : Fraction of the rudder angle to add to the angle of attack.
        Must be between 0 and 1.
        Normally, use 1. for a rudder, 0. otherwise.
    speed_multiplier : Multiplier of the free flow speed.
        Must be between 0.6 and 1.4.
        A current practice is, for example, to use a 0.9 speed multiplier for the rudder.
    rho_water : Water density [kg/m**3], must be >= 0

    Returns a Force object, representing the sideforce and induced drag [N] on the hull.

    References
    ----------
    Aero-Hydrodynamics and the performance of sailing yachts - Fabio Fossati - p143/144

    """
    if projected_area < 0:
        raise ValueError("projected_area should be positive or zero")
    if span <= 0:
        raise ValueError("span should be strictly positive")
    if not 0.5 <= aspect_ratio_multiplier <= 2.5:
        raise ValueError("unrealistic aspect_ratio_multiplier")
    if not 0 <= rudder_angle_influence <= 1.:
        raise ValueError("unrealistic rudder_angle_influence")
    if not 0.6 <= speed_multiplier <= 1.4:
        raise ValueError("unrealistic speed_multiplier")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    if filename_2d_section is None and add_viscous_drag is True:
        raise ValueError("Cannot compute viscous Cd without a 2d foil dat file")

    geometric_aspect_ratio = span**2 / projected_area
    effective_aspect_ratio = geometric_aspect_ratio * aspect_ratio_multiplier
    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.

    if use_cos_heel_angle is True:
        aoa = (leeway_angle + rudder_angle_influence * rudder_angle) * cos(radians(heel_angle))
    else:
        aoa = (leeway_angle + rudder_angle_influence * rudder_angle)

    if filename_2d_section is None:
        slope = 0.105
    else:
        raise NotImplementedError

    cl_alpha_3d = slope / (1. + (1.92 / effective_aspect_ratio))

    Cl = cl_alpha_3d * aoa

    lift = Cl * 0.5 * rho_water * projected_area * (boatspeed * speed_multiplier)**2

    # induced drag coefficient
    oswald_efficiency_factor = 0.9
    Cdi = Cl**2 / (pi * effective_aspect_ratio * oswald_efficiency_factor)

    induced_drag = 0.5 * Cdi * rho_water * projected_area * (boatspeed * speed_multiplier)**2

    if add_viscous_drag is True:
        # 2d section dat file viscous drag
        raise NotImplementedError
    else:
        # 3d induced drag only
        drag = induced_drag

    px, py, pz = shift_coe_around_x(coe, heel_angle)

    return Force(-drag * boatspeed_sign,
                 -lift * cos(radians(heel_angle)) * boatspeed_sign,
                 lift * sin(radians(heel_angle)) * boatspeed_sign,
                 px, py, pz)


def beta_zero(beam_to_draft_ratio: float, heel_angle_degrees: float) -> float:
    """zero lift drift angle. Is always positive !!"""
    return (0.405 * beam_to_draft_ratio * radians(heel_angle_degrees)) ** 2


def keel_downwash_on_rudder(keel_Cl: float,
                            keel_effective_aspect_ratio: float,
                            heel_angle_degrees: float) -> float:
    """Is always positive !!"""
    # if keel_Cl / keel_effective_aspect_ratio < 0.:
    #   return 0.
    # else:
    #   return degrees((0.136 + 0.001 * abs(heel_angle_degrees) / 15.)
    #                            * sqrt( keel_Cl / keel_effective_aspect_ratio))
    # if keel_effective_aspect_ratio == 0:
    #     raise ValueError
    # else:
    return degrees((0.136 + 0.001 * abs(heel_angle_degrees) / 15.) * sqrt(abs(keel_Cl) / keel_effective_aspect_ratio))


def lift_curve_slope(effective_aspect_ratio: float, sweep_back_angle: float) -> float:
    """Lift curve slope for an attack angle in degrees"""
    lift_curve_slope_for_radians = (5.7 * effective_aspect_ratio
                                    / (1.8 + (cos(radians(sweep_back_angle))
                                              * sqrt(effective_aspect_ratio ** 2
                                                     / (cos(radians(sweep_back_angle)) ** 4) + 4.))))
    lift_curve_slope_for_degrees = lift_curve_slope_for_radians * pi / 180.
    return lift_curve_slope_for_degrees


def sideforce_production_keuning_verwerft(boatspeed: float,
                                          heel_angle: float,
                                          leeway_angle: float,
                                          rudder_angle: float,
                                          keel_projected_area: float,
                                          keel_span: float,
                                          rudder_projected_area: float,
                                          rudder_span: float,
                                          lwl: float,
                                          bwl: float,
                                          Tc: float,
                                          hull_coe: Optional[Tuple[float, float, float]] = None,
                                          keel_coe: Optional[Tuple[float, float, float]] = None,
                                          rudder_coe: Optional[Tuple[float, float, float]] = None,
                                          keel_sweep_back_angle: float = 0.,
                                          rudder_sweep_back_angle: float = 0.,
                                          rho_water: float = RHO_SEA_WATER_20C) -> Tuple[Force, Force, Force]:
    """Sideforce the Keuning Verwerft way.

    Calculate the lift (i.e. sideforce) and drag due to lift [N] according
    to the 2009 publication by Keuning and Verwerft.
    Viscous effects are not accounted for;
    other components are available to represent the viscous drag of the hull and its appendages.

    boatspeed : boat speed [m/s]
    heel_angle : The heel_angle [degrees], between -90 and 90.
    leeway_angle : The leeway angle [degrees]
    rudder_angle : The rudder angle [degrees]
    keel_projected_area : Projected area of the keel [m**2], must be >= 0
    keel_span : Span of the keel [m] , must be >= 0.
    rudder_projected_area : Projected area of the rudder [m**2], must be >= 0
    rudder_span : Span of the rudder [m] , must be >= 0.
    lwl : Length at Waterline [m] , must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    Tc : Canoe body draft [m], must be >= 0
    hull_coe : tuple of 3 floats, optional (default [lwl * 3/4, 0., -Tc / 3])
        Represents the x, y, z coordinates of the centre of effort.
    keel_coe : tuple of 3 floats
               optional (default [0., 0., -(Tc + keel_span) * 0.43])
        NOT PROVIDING A COE MAKE SOLVING BEYOND 3 DOF IMPOSSIBLE.
        Represents the x, y, z coordinates of the centre of effort.
        The centre of effort can be estimated by projecting the projected
        surface barycentre to the 1/4 chord line.
    rudder_coe : tuple of 3 floats
                 optional (default [0., 0., - rudder_span * 0.43])
        NOT PROVIDING A COE MAKE SOLVING BEYOND 3 DOF IMPOSSIBLE.
        Represents the x, y, z coordinates of the centre of effort.
        The centre of effort can be estimated by projecting the
        projected surface barycentre to the 1/4 chord line.
    keel_sweep_back_angle : The sweep back angle [deg.] of the 1/4 chord line of the keel.
        Must be between -45 and 45
    rudder_sweep_back_angle : The sweep back angle [deg.] of the 1/4 chord line of the rudder.
        Must be between -45 and 45
    rho_water : Water density [kg/m**3], must be >= 0

    Returns a list of Force objects, representing the sideforce and induced drag [N]
    of the keel, rudder and hull in this order

    References
    ----------
    A new Method for the Prediction of the Side Force on Keel and Rudder
    of a Sailing Yacht based on the Results of the Delft Systematic Yacht
    Hull Series by J. A. Keuning and B. Verwerft (March 2009)
    File : CSYS19Keuning.pdf (web name)

    """
    if keel_projected_area < 0:
        raise ValueError("keel_projected_area should be positive or zero")
    if keel_span <= 0:
        raise ValueError("keel_span should be strictly positive")
    if rudder_projected_area < 0:
        raise ValueError("rudder_projected_area should be positive or zero")
    if rudder_span <= 0:
        raise ValueError("rudder_span should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if not -45. <= keel_sweep_back_angle <= 45.:
        raise ValueError("unrealistic keel_sweep_back_angle")
    if not -45. <= rudder_sweep_back_angle <= 45.:
        raise ValueError("unrealistic rudder_sweep_back_angle")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    keel_geometric_aspect_ratio = keel_span ** 2 / keel_projected_area
    rudder_geometric_aspect_ratio = rudder_span ** 2 / rudder_projected_area

    keel_effective_aspect_ratio = 2. * keel_geometric_aspect_ratio
    rudder_effective_aspect_ratio = 2 * rudder_geometric_aspect_ratio

    # lift_curve_slope is the lift per unit angle of attack (Whicker and Fehlner)
    keel_lift_curve_slope = lift_curve_slope(keel_effective_aspect_ratio,
                                             keel_sweep_back_angle)
    rudder_lift_curve_slope = lift_curve_slope(rudder_effective_aspect_ratio,
                                               rudder_sweep_back_angle)

    # hull influence coefficient
    c_hull = (1.80 * Tc / keel_span) + 1

    if hull_coe is None:
        # keel_coe = [0., 0., -Tc - keel_span / 2.]
        hull_coe = (lwl * 3. / 4., 0., -Tc / 3.)

    if keel_coe is None:
        # keel_coe = [0., 0., -Tc - keel_span / 2.]
        keel_coe = (0., 0., -(Tc + keel_span) * 0.43)

    if rudder_coe is None:
        # rudder_coe = [0., 0., - rudder_span / 2.]
        rudder_coe = (0., 0., - rudder_span * 0.43)

    # Signs
    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.
    leeway_angle_sign = leeway_angle / abs(leeway_angle) if leeway_angle != 0. else 0.
    heel_angle_sign = heel_angle / abs(heel_angle) if heel_angle != 0. else 0.

    # make sure no solution can be found for leeway > 15.
    leeway_angle = leeway_angle if abs(leeway_angle) < 15. else 0.

    # heel influence coefficient (c_heel = 1 @ 0° heel, c_heel = 0.4 @ 90° heel)
    c_heel = 1 - 0.382 * radians(abs(heel_angle))

    angle_of_attack_on_keel = ((leeway_angle -
                                (beta_zero(bwl / Tc, heel_angle) *
                                 heel_angle_sign)) * cos(radians(heel_angle)))

    keel_Cl = keel_lift_curve_slope * angle_of_attack_on_keel
    # print('keel : %s' % str(keel_Cl))

    # downwash angle of the keel on the rudder
    # keel_downwash_on_rudder = degrees((0.136 + 0.001 * heel_angle / 15.)
    #                           * sqrt( keel_Cl / keel_effective_aspect_ratio))
    keel_downwash = keel_downwash_on_rudder(keel_Cl,
                                            keel_effective_aspect_ratio,
                                            heel_angle)

    angle_of_attack_on_rudder = ((leeway_angle -
                                  (beta_zero(bwl / Tc, heel_angle) *
                                   heel_angle_sign)
                                  - (keel_downwash * leeway_angle_sign) +
                                  rudder_angle)
                                 * cos(radians(heel_angle)))
    rudder_Cl = rudder_lift_curve_slope * angle_of_attack_on_rudder
    # print('rudder : %s' % str(rudder_Cl))

    lift_keel = (c_hull * c_heel * 0.5 * rho_water * keel_projected_area *
                 boatspeed ** 2 * keel_Cl)
    # * leeway_angle_sign
    lift_rudder = (c_hull * c_heel * 0.5 * rho_water * rudder_projected_area *
                   (0.9 * boatspeed) ** 2 * rudder_Cl)
    # * leeway_angle_sign

    # 0.9 is an approximation of wing span efficiency
    # 15.  Perkins, C.D. and Ciale, R.E.,
    #                             "Airplane  Performance Stability and Control,"
    # JohnWiley and Sons, Inc.  (1949).
    oswald_efficiency_factor = 0.9
    induced_drag_keel = (lift_keel ** 2
                         /
                         (0.5 * rho_water * keel_projected_area *
                          boatspeed ** 2 * pi * oswald_efficiency_factor *
                          keel_effective_aspect_ratio)) if boatspeed != 0. else 0.

    induced_drag_rudder = (lift_rudder ** 2
                           /
                           (0.5 * rho_water * rudder_projected_area *
                            (0.9 * boatspeed) ** 2 * pi
                            * oswald_efficiency_factor *
                            rudder_effective_aspect_ratio)) if boatspeed != 0. else 0.

    lift_hull = (lift_keel + lift_rudder) - ((lift_keel + lift_rudder) / c_hull)
    # induced_drag = induced_drag_keel + induced_drag_rudder

    px_k, py_k, pz_k = shift_coe_around_x(keel_coe, heel_angle)

    keel_force = Force(-induced_drag_keel * boatspeed_sign,
                       -lift_keel * boatspeed_sign / c_hull,
                       lift_keel * tan(radians(heel_angle)) * boatspeed_sign,
                       px_k, py_k, pz_k)

    px_r, py_r, pz_r = shift_coe_around_x(rudder_coe, heel_angle)

    rudder_force = Force(-induced_drag_rudder * boatspeed_sign,
                         -lift_rudder * boatspeed_sign / c_hull,
                         (lift_rudder * tan(radians(heel_angle)) * boatspeed_sign),
                         px_r, py_r, pz_r)

    px_h, py_h, pz_h = shift_coe_around_x(hull_coe, heel_angle)

    hull_force = Force(0., -lift_hull * boatspeed_sign, 0., px_h, py_h, pz_h)

    return keel_force, rudder_force, hull_force
