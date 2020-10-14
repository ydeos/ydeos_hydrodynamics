# coding: utf-8

r"""Viscous resistance estimates."""

from typing import Tuple, Callable, Optional
from ydeos_hydrodynamics.constants import RHO_SEA_WATER_20C, KINEMATIC_VISCOSITY_WATER_20C
from ydeos_hydrodynamics.friction_lines import ittc57, hughes, prandtl
from ydeos_hydrodynamics.form_factors import hoerner, hoerner_streamlined_body
from ydeos_hydrodynamics.geometry import shift_coe_around_x
from ydeos_hydrodynamics.force import Force
from ydeos_hydrodynamics.water import RHO_WATER_MIN, RHO_WATER_MAX


def hull_viscous(boatspeed: float,
                 lwl: float,
                 Sc: float,
                 coe: Optional[Tuple[float, float, float]] = None,
                 lwl_multiplier: float = 0.7,
                 transition_reynolds_number: float = 5e5,
                 form_factor: float = 1.12,
                 rho_water: float = RHO_SEA_WATER_20C,
                 kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
                 turbulent_friction_line: Callable = ittc57) -> Force:
    """Calculate the friction resistance [N] on the hull.

    boatspeed : boat speed [m/s]
    lwl : Length at Waterline [m] , must be > 0.
    Sc : WSAc - Canoe body Wetted Surface Area [m**2], must be >= 0
    coe : tuple of 3 floats, optional (default [lwl / 2, 0, -Tc/3]).
        Represents the x, y, z coordinates of the centre of effort.
        For viscous resistance, the surfacic barycentre of the underwater
            hull is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    lwl_multiplier : Multiplier of the waterline length
        (used for Reynolds number calculation).
        The ITTC 57 practice is to use 0.7. Must be between 0.7 and 1.0.
    transition_reynolds_number : The Reynolds number at which the flow
        transitions from laminar to turbulent.
        Must be >=0. A value of 0 denotes a fully turbulent flow.
    form_factor : form factor to account for pressure drag
        (default 1.12 (embedded ITTC57 form factor))
        Must be between 0.9 and 1.6
    rho_water : Water density [kg/m**3], must be >= 0.
    kinematic_viscosity : defaults to water at 20°C : 1,004*10E-6 [m**2/s]
    turbulent_friction_line : the function used to compute the turbulent Cf
        (e.g. ittc57, hughes), used in prandtl

    Returns a Force object, representing the friction resistance [N] on the hull.

    References
    ----------
    http://radiosailingtechnology.com/index.php/hulls/drag-measurements-on-an-international-one-metre-yacht

    """
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if Sc <= 0:
        raise ValueError("Sc should be strictly positive")
    if Sc >= lwl * lwl / 2:
        raise ValueError("Sc is unrealistic")
    if not 0.7 <= lwl_multiplier <= 1:
        raise ValueError("unrealistic lwl_multiplier")
    if transition_reynolds_number < 0:
        raise ValueError("transition_reynolds_number should be positive or zero")
    if not 0.9 <= form_factor <= 1.6:
        raise ValueError("unrealistic form_factor")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")
    if kinematic_viscosity <= 0:
        raise ValueError("kinematic viscosity should be strictly positive")
    if turbulent_friction_line not in [ittc57, hughes]:
        raise ValueError("Unknown turbulent friction line")

    if coe is None:
        coe = (lwl / 2., 0., 0.)

    boatspeed_sign = 1. if boatspeed == 0. else boatspeed / abs(boatspeed)

    # TODO: where does the division of the form factor by 1.12 come from?
    Cf = (form_factor / 1.12) * prandtl(boatspeed,
                                        lwl_multiplier * lwl,
                                        transition_reynolds_number,
                                        kinematic_viscosity,
                                        turbulent_friction_line)

    friction_resistance = 0.5 * rho_water * boatspeed ** 2 * Sc * Cf

    return Force(-friction_resistance * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])


def appendage_viscous(boatspeed: float,
                      heel_angle: float,
                      coe: Tuple[float, float, float],
                      average_chord: float,
                      thickness_to_chord: float,
                      wetted_area: float,
                      transition_reynolds_number: float = 5e5,
                      rho_water: float = RHO_SEA_WATER_20C,
                      kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
                      turbulent_friction_line: Callable = ittc57) -> Force:
    """Calculate the viscous resistance [N] on an appendage.

    (For example, on a keel fin, a rudder, a centreboard ....)
    The centre of effort position is approximated to rotate
    around the X axis with heel.

    boatspeed : boat speed [m/s]
    heel_angle : The heel_angle [degrees], between -90 and 90.
    coe :  tuple of 3 floats
        Represents the x, y, z coordinates of the centre of effort.
        For viscous resistance, the surfacic barycentre of the appendage
        is an acceptable estimate.
    average_chord : Average chord of the appendage [m] , must be > 0.
    thickness_to_chord : Thickness to chord ratio of keel
        (e.g.0.1 for NACA0010). Must be between 0 and 0.4
    wetted_area : WSA - Wetted Surface Area of the appendage [m**2],
                  must be >= 0
    transition_reynolds_number : The Reynolds number at which the flow
        transitions from laminar to turbulent.
        A value of 0 denotes a fully turbulent flow.
        Optional (default 500 000 (5e5)). Must be >=0.
    rho_water : Water density [kg/m**3], must be >= 0 (default RHO_SEA_WATER)
    kinematic_viscosity : defaults to water at 20°C : 1,004*10E-6 [m**2/s]
    turbulent_friction_line : the function used to compute the turbulent Cf
        (e.g. ittc57, hughes), used in prandtl

    Returns a Force object, representing the viscous resistance [N]
    on an appendage.

    References
    ----------
    www.boatdesign.net/forums/attachments/hydrodynamics-aerodynamics/
    79702d1363296908-keel-rudder-viscous-drag-correlation-line.pdf

    """
    if average_chord <= 0:
        raise ValueError("average_chord should be strictly positive")
    if not 0. <= thickness_to_chord < 0.4:
        raise ValueError("unrealistic thickness_to_chord")
    if wetted_area <= 0:
        raise ValueError("wetted_area should be strictly positive")
    if transition_reynolds_number < 0:
        raise ValueError("transition_reynolds_number should be positive or zero")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between "
                         f"{RHO_WATER_MIN} and {RHO_WATER_MAX}")
    if kinematic_viscosity <= 0:
        raise ValueError("kinematic viscosity should be strictly positive")
    if turbulent_friction_line not in [ittc57, hughes]:
        raise ValueError("Unknown turbulent friction line")

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.

    cf = prandtl(boatspeed, average_chord, transition_reynolds_number,
                 kinematic_viscosity, turbulent_friction_line)

    viscous_resistance = hoerner(thickness_to_chord) * 0.5 * rho_water * boatspeed ** 2 * wetted_area * cf
    px, py, pz = shift_coe_around_x(coe, heel_angle)

    return Force(-viscous_resistance * boatspeed_sign, 0., 0., px, py, pz)


# Ballast (streamlined) resistance functions


def ballast_viscous(boatspeed: float,
                    heel_angle: float,
                    coe: Tuple[float, float, float],
                    bulb_length: float,
                    bulb_diameter: float,
                    bulb_wetted_area: float,
                    transition_reynolds_number: float = 5e5,
                    rho_water: float = RHO_SEA_WATER_20C,
                    kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
                    turbulent_friction_line: Callable = ittc57) -> Force:
    """Calculates the viscous resistance [N] on a bulb.

    The centre of effort position is approximated to rotate
    around the X axis with heel.

    boatspeed : boat speed [m/s]
    heel_angle : The heel_angle [degrees], between -90 and 90.
    coe : tuple of 3 floats
        Represents the x, y, z coordinates of the centre of effort.
        For viscous resistance, the surfacic barycentre of the bulb
        is an acceptable estimate.
    bulb_length : length of the bulb [m] , must be > 0.
    bulb_diameter : Maximum diameter [m] of the bulb,
                   if the cross section is circular
                   Otherwise, use : sqrt (4 * Amax / pi) where Amax is the
                   maximum cross section area
                   bulb_diameter / bulb_length must be between 0 and 0.4
    bulb_wetted_area : Total Bulb Wetted Surface Area [m**2], must be >= 0
    transition_reynolds_number : The Reynolds number at which the flow
        transitions from laminar to turbulent.
        (default 500 000 (5e5)). Must be >=0.
        A value of 0 denotes a fully turbulent flow.
    rho_water : Water density [kg/m**3], must be >= 0 (default RHO_SEA_WATER)
    kinematic_viscosity : defaults to water at 20°C : 1,004*10E-6 [m**2/s]
    turbulent_friction_line : the function used to compute the turbulent Cf
        (e.g.    ittc57, hughes), used in prandtl

    Returns a Force object, representing the viscous resistance [N]
    on an appendage.

    References
    ----------
    http://www.boatdesign.net/forums/boat-design/keel-bulb-form-factors-4569.html
    http://books.google.fr/books?id=XtU4HVnWeZIC&pg=PA702#v=onepage&q&f=false

    """
    if bulb_length <= 0:
        raise ValueError("bulb_length should be strictly positive")
    if bulb_diameter <= 0:
        raise ValueError("bulb_diameter should be strictly positive")
    if bulb_diameter / bulb_length >= 0.4:
        raise ValueError("The bulb_diameter / bulb_length ratio is unrealistic")
    if bulb_wetted_area <= 0:
        raise ValueError("bulb_wetted_area should be strictly positive")
    if transition_reynolds_number < 0:
        raise ValueError("transition_reynolds_number should be positive or zero")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between "
                         f"{RHO_WATER_MIN} and {RHO_WATER_MAX}")
    if kinematic_viscosity <= 0:
        raise ValueError("kinematic viscosity should be strictly positive")
    if turbulent_friction_line not in [ittc57, hughes]:
        raise ValueError("Unknown turbulent friction line")

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.

    c_friction = prandtl(boatspeed, bulb_length, transition_reynolds_number,
                         kinematic_viscosity, turbulent_friction_line)

    form_factor = hoerner_streamlined_body(bulb_diameter / bulb_length)

    viscous_resistance = form_factor * 0.5 * rho_water * boatspeed ** 2 * bulb_wetted_area * c_friction
    px, py, pz = shift_coe_around_x(coe, heel_angle)

    return Force(-viscous_resistance * boatspeed_sign, 0., 0., px, py, pz)
