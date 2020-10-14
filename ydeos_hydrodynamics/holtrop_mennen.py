# coding: utf-8

r"""Holtrop-Mennen resistance estimates."""

from typing import Tuple
from math import sqrt, exp, pi, cos
from ydeos_hydrodynamics.constants import KINEMATIC_VISCOSITY_WATER_20C, \
    RHO_SEA_WATER_20C, GRAVITY_STANDARD
from ydeos_hydrodynamics.friction_lines import ittc57
from ydeos_hydrodynamics.froude import froude_number


def holtrop_wetted_area(lwl: float,
                        B: float,
                        T: float,
                        Cm: float,
                        Cb: float,
                        Cwp: float,
                        Abt: float) -> float:
    r"""Wetted area estimate from hull parameters.

    Cb : block coefficient
    Cwp : coefficient of waterplane
    Abt : transverse bulb area

    """
    return lwl * (2 * T + B) * sqrt(Cm) * (0.453 + 0.4425 * Cb - 0.2862 * Cm - 0.003467 * B / T + 0.3696 * Cwp) + 2.38 * Abt / Cb


def _lcb_fraction_to_lcb_holtrop(lcb_fraction: float) -> float:
    r"""0.5x format conversion to holtrop percentage."""
    return (0.5 - lcb_fraction) * 100


def _holtrop_Lr(lwl: float, Cp: float, lcb_fraction: float) -> float:
    r"""Run length."""
    return lwl * (1 - Cp + 0.06 * Cp * _lcb_fraction_to_lcb_holtrop(lcb_fraction) / (4 * Cp - 1))


def _c3func(Abt, B, T, Tf, h_B):
    r"""Holtrop c3 coefficient."""
    return 0.56 * Abt**1.5 / (B * T * (0.31 * sqrt(Abt) + Tf - h_B))


def _c2func(c3_):
    r"""Holtrop c2 coefficient.

    c2 is a parameter that accounts for
    the reduction of the wave resistance due to the action of a bulbous bow

    """
    return exp(-1.89 * sqrt(c3_))


def holtrop_friction(v: float,
                     lwl: float,
                     B: float,
                     T: float,
                     Cp: float,
                     lcb_fraction: float,
                     afterbody_shape: str,  # U N or V
                     S: float,
                     rho_water: float = RHO_SEA_WATER_20C,
                     kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> Tuple[float, ...]:
    r"""Holtrop friction drag and form factor."""
    lcb = _lcb_fraction_to_lcb_holtrop(lcb_fraction)
    Lr = _holtrop_Lr(lwl, Cp, lcb_fraction)
    c12 = (T / lwl)**0.2228446 if T / lwl > 0.05 else 0.479948 if T / lwl < 0.02 else 48.2 * (T/lwl - 0.02)**2.078 + 0.479948
    c13 = 1 + 0.003 * {"V": -10, "N": 0, "U": 10}[afterbody_shape]
    k1 = (c13 * (0.93 + c12 * (B / Lr)**0.92497 * (0.95 - Cp)**(-0.521448) * (1 - Cp + 0.0225 * lcb)**0.6906)) - 1
    Cf = ittc57(v, lwl, kinematic_viscosity)
    return 0.5 * rho_water * S * v**2 * Cf, k1, c12, c13, Cf


def holtrop_appendages(v: float,
                       S_app: Tuple[float, ...] = (),
                       d_app: Tuple[float, ...] = (),
                       k2_app: Tuple[float, ...] = (),
                       rho_water: float = RHO_SEA_WATER_20C,
                       kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C) -> float:
    r"""Holtrop resistance of appendages.

    S_app : tuple of appendage surfaces
    d_app : tuple of appendage typical dimension
    k2_app: tuple of appendage k2 values

    Approximate 1+k2 values
    -----------------------
    rudder behind skeg          1.5 - 2.0
    rudder behind stern         1.3 - 1.5
    twin-screw balance rudder   2.8
    shaft brackets              3.0
    skeg                        1.5 - 2.0
    strut bossings              3.0
    hull bossings               2.0
    shafts                      2.0 - 4.0
    stabilizer fins             2.8
    dome                        2.7
    bilge keels                 1.4

    """
    if not len(S_app) == len(d_app) == len(k2_app):
        raise ValueError("Different lengths for S_app, d_app and k2_app")

    r_appendages = sum(
        [(lambda d, S, k2: 0.5 * rho_water * S * v ** 2 * (1 + k2) * ittc57(v, d,
                                                                            kinematic_viscosity))(
            t[1], t[0], t[2])
         for t in zip(S_app, d_app, k2_app)])
    return r_appendages


def holtrop_bow_thruster_opening(v: float,
                                 d_bto: float,
                                 Cbto: float = 0.0075,
                                 rho_water: float = RHO_SEA_WATER_20C) -> float:
    r"""Bow thruster opening resistance.

    not included in the global Holtrop-Mennen method.

    d_bto : bow thruster opening diameter
    Cbto: coefficient from 0.003 to 0.012

    """
    return rho_water * v ** 2 * pi * d_bto ** 2 * Cbto


def holtrop_wave(v: float,
                 Vc: float,
                 lwl: float,
                 B: float,
                 T: float,
                 Cp: float,
                 lcb_fraction: float,
                 At: float,  # area transom
                 Cm: float,
                 Cwp: float,  # waterplane area coefficient
                 Abt: float,  # Transverse bulb area
                 h_B: float,
                 Tf: float,  # forward draught of the ship
                 rho_water: float = RHO_SEA_WATER_20C,
                 gravity: float = GRAVITY_STANDARD) -> Tuple[float, ...]:
    r"""Holtrop wave-making and wave-breaking resistance.

    At : immersed part of the transverse area of the transom  at zero speed

    h_B : position of the centre of the transverse area Abt above the keel line
    Tf : forward draught of the ship

    """
    Fn = froude_number(v, lwl, gravity)
    Lr = _holtrop_Lr(lwl, Cp, lcb_fraction)
    lcb = _lcb_fraction_to_lcb_holtrop(lcb_fraction)

    # The half angle of entrance i_e is the angle of the waterline
    # at bow in degrees with reference to the centreplane
    i_e = 1 + 89 * exp(
        -(lwl / B) ** 0.80856 * (1 - Cwp) ** 0.30484 * (1 - Cp - 0.0225 * lcb) ** 0.6367 * (
                    Lr / B) ** 0.34574 * (100 * Vc / lwl ** 3) ** 0.16302)

    c7 = 0.229577 * (B / lwl) ** 0.33333 if B / lwl < 0.11 else 0.5 - 0.0625 * (
                lwl / B) if B / lwl > 0.25 else B / lwl
    c1 = 2223105 * c7 ** 3.78613 * (T / B) ** 1.07961 * (90 - i_e) ** (-1.37565)
    c3 = _c3func(Abt, B, T, Tf, h_B)
    c2 = _c2func(c3)

    # c5 expresses the influence of a transom stern on the resistance
    c5 = 1 - 0.8 * At / (B * T * Cm)

    c16 = 8.07981 * Cp - 13.8673 * Cp ** 2 + 6.984388 * Cp ** 3 if Cp < 0.8 else 1.73014 - 0.7067 * Cp
    m1 = 0.0140407 * lwl / T - 1.75254 * Vc ** (1. / 3) / lwl - 4.79323 * B / lwl - c16
    lambda_ = 1.446 * Cp - 0.03 * (lwl / B) if lwl / B < 12 else 1.446 * Cp - 0.36
    c15 = -1.69385 if lwl ** 3 / Vc < 512 else 0. if lwl ** 3 / Vc > 1727 else -1.69385 + (
                lwl / Vc ** (1. / 3) - 8.0) / 2.36
    m2 = c15 * Cp ** 2 * exp(-0.1 * Fn ** (-2))
    d = -0.9
    # TODO: cos ou cos(radians?
    return c1 * c2 * c5 * Vc * rho_water * gravity * exp(
        m1 * Fn ** d + m2 * cos(lambda_ * Fn ** (-2))), c1, c5, c7, c15, m1, m2, lambda_


def holtrop_bulbous(v: float,
                    Abt: float,
                    h_B: float,
                    Tf: float,
                    rho_water: float = RHO_SEA_WATER_20C,
                    gravity: float = GRAVITY_STANDARD) -> Tuple[float, ...]:
    r"""Holtrop additional resistance of bulbous bow near the water surface.

    Parameters
    ----------
    Abt : Transverse bulb area
    h_B : h_B is the position of the centre of the transverse area Abt above the keel line
    Tf : forward draught of the ship

    """
    Pb = 0.56 * sqrt(Abt) / (Tf - 1.5 * h_B)
    Fni = v / sqrt(gravity * (Tf - h_B - 0.25 * sqrt(Abt)) + 0.15 * v**2)

    return 0.11 * exp(-3 * Pb**(-2)) * Fni**3 * Abt**1.5 * rho_water * gravity / (1 + Fni**2), Pb, Fni


def holtrop_transom(v: float,
                    B: float,
                    At: float,  # area transom
                    Cwp: float,  # waterplane area coefficient
                    rho_water: float = RHO_SEA_WATER_20C,
                    gravity: float = GRAVITY_STANDARD) -> Tuple[float, ...]:
    r"""Holtrop additional pressure resistance of immersed transom stern."""
    FnT = v / sqrt(2 * gravity * At / (B + B * Cwp))
    c6 = 0.2 * (1 - 0.2 * FnT) if FnT < 0.5 else 0
    return 0.5 * rho_water * v**2 * At * c6, c6, FnT


def holtrop_a(v: float,
              lwl: float,
              B: float,
              T: float,
              S: float,
              Cb: float,
              Abt: float,
              h_B: float,
              Tf: float,
              rho_water: float = RHO_SEA_WATER_20C) -> Tuple[float, ...]:
    r"""Holtrop model-ship correlation resistance.

    Primarily describes the effect of the hull roughness
    and the still-air resistance

    Parameters
    ----------
    Cb : block coefficient
    Abt : Transverse bulb area
    h_B : h_B is the position of the centre of the transverse area Abt
          above the keel line
    Tf : forward draught of the ship

    """
    c3 = _c3func(Abt, B, T, Tf, h_B)
    c2 = _c2func(c3)
    c4 = Tf / lwl if Tf/lwl <= 0.04 else 0.04
    Ca = 0.006 * (lwl + 100)**(-0.16) - 0.00205 + 0.003 * sqrt(lwl/7.5) * Cb**4 * c2 * (0.04 - c4)
    return 0.5 * rho_water * v**2 * S * Ca, c4, Ca


def holtrop_mennen(v: float,
                   Vc: float,
                   lwl: float,
                   B: float,
                   T: float,
                   Cp: float,
                   lcb_fraction: float,
                   afterbody_shape: str,
                   S: float,
                   At: float,  # area transom (immersed)
                   Cm: float,
                   Cb: float,  # block
                   Cwp: float,  # waterplane area coefficient
                   Abt: float,  # Transverse bulb area
                   h_B: float,  # h_B is the position of the centre of the transverse area Abt above the keel line
                   Tf: float,  # forward draught of the ship
                   S_app: Tuple[float, ...] = (),
                   d_app: Tuple[float, ...] = (),
                   k2_app: Tuple[float, ...] = (),
                   rho_water: float = RHO_SEA_WATER_20C,
                   kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
                   gravity: float = GRAVITY_STANDARD) -> Tuple[Tuple[float, float, float, float, float, float, float], float]:
    r"""Holtrop Mennen model.

    Holtrop-Mennen power prediction of high block ships with low L/B ratios
    and of slender naval ships with a complex appendage arrangement
    and immersed transom sterns

    The output is not wrapped in a Force since it will probably
    never be used in a VPP.

    h_B : position of the centre of the transverse area Abt above the keel line
    Tf : forward draught of the ship

    Limitations
    -----------
    Ship type	            No. Froude máx.	    Cp min-max	L/B min-max	    B/T min-max

    Tankers, Bulk carriers	0,24	            0,73-0,85	5,1-7,1	        2,4-3,2
    Trawler, Coastal, Tug	0,38	            0.55-0.65	3.9-6.3	        2.1-3.0
    Containers  	        0.45	            0.55-0.67	6.0-9.5	        3.0-4.0
    Cargo   	            0.30	            0.56-0.75	5.3-8.0	        2.4-4.0
    Ro-ro, ferries          0.35	            0.55-0.67	5.3-8.0	        3.2-4.0

    """
    r_friction, k1, c12, c13, Cf = holtrop_friction(v,
                                                    lwl,
                                                    B,
                                                    T,
                                                    Cp,
                                                    lcb_fraction,
                                                    afterbody_shape,
                                                    S,
                                                    rho_water,
                                                    kinematic_viscosity)
    r_appendages = holtrop_appendages(v,
                                      S_app,
                                      d_app,
                                      k2_app,
                                      rho_water,
                                      kinematic_viscosity)
    r_wave, c1, c5, c7, c15, m1, m2, lambda_ = holtrop_wave(v,
                                                            Vc,
                                                            lwl,
                                                            B,
                                                            T,
                                                            Cp,
                                                            lcb_fraction,
                                                            At,
                                                            Cm,
                                                            Cwp,
                                                            Abt,
                                                            h_B,
                                                            Tf,
                                                            rho_water,
                                                            gravity)
    r_bulbous, Pb, Fni = holtrop_bulbous(v, Abt, h_B, Tf, rho_water, gravity)
    r_transom, c6, FnT = holtrop_transom(v, B, At, Cwp, rho_water, gravity)
    r_a, c4, Ca = holtrop_a(v, lwl, B, T, S, Cb, Abt, h_B, Tf, rho_water)

    total = (1 + k1) * r_friction + r_appendages + r_wave + r_bulbous + r_transom + r_a

    return (r_friction, (1 + k1), r_appendages, r_wave, r_bulbous, r_transom, r_a), total
