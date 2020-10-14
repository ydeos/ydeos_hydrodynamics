# coding: utf-8

r"""Reynolds number computations."""

from ydeos_hydrodynamics.constants import KINEMATIC_VISCOSITY_WATER_20C


def reynolds_number(speed: float,
                    dimension: float,
                    kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> float:
    """Reynolds number.

    Reynolds number from flow speed [m/s], dimension [m]
    and kinematic viscosity [m2/s].

    kinematic_viscosity defaults to kinematic viscosity of water at 20°C :
    1,004*10E-6 (m2/s) at 20°C [m**2/s]

    """
    if dimension < 0:
        raise ValueError("dimension should be positive or zero")
    if kinematic_viscosity <= 0:
        raise ValueError("kinematic_viscosity should be strictly positive")
    return (abs(speed) * dimension) / kinematic_viscosity


def transition_x(speed: float,
                 transition_reynolds_number: float,
                 kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> float:
    """Transition position.

    Distance [m] from the leading edge at which the local Reynolds number
    is equal to the transition_reynolds_number.

    speed [m/s]
    transition_reynolds_number is the Reynolds number at which the transition
    from laminar to turbulent is expected to occur
    kinematic_viscosity defaults to  water at 20°C : 1,004*10E-6 [m**2/s]

    """
    if speed == 0.:
        # return sys.float_info.max
        return float("inf")  # infinity
    return transition_reynolds_number * kinematic_viscosity / abs(speed)
