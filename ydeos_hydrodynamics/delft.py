# coding: utf-8

r"""Delft resistance estimates"""

from typing import Tuple, Union
from math import radians, sqrt
from functools import partial
import numpy as np
from scipy.interpolate import interp1d
from ydeos_hydrodynamics.domain import check_domain
from ydeos_hydrodynamics.constants import RHO_SEA_WATER_20C, GRAVITY_STANDARD
from ydeos_hydrodynamics.force import Force
from ydeos_hydrodynamics.water import RHO_WATER_MIN, RHO_WATER_MAX
from ydeos_hydrodynamics.memoize import memoize
from ydeos_hydrodynamics.froude import froude_number

# Geometry checks against the Delft Systematic Yacht Hull Series limitations
# A warning is issued if the limitations are exceeded in the check* procedures

DSYHS_2008 = "DSYHS 2008"
DSYHS_S4 = "DSYHS S4"
DSYHS_S3 = "DSYHS S3"

DSYHS_S3_BEAM_TO_DRAFT_MIN = DSYHS_S4_BEAM_TO_DRAFT_MIN = DSYHS_2008_BEAM_TO_DRAFT_MIN = 2.46
DSYHS_S4_BEAM_TO_DRAFT_MAX = DSYHS_2008_BEAM_TO_DRAFT_MAX = 19.38
DSYHS_S3_BEAM_TO_DRAFT_MAX = 19.32

DSYHS_S4_WATERPLANE_TO_VOLUME_MIN = 3.78
DSYHS_S4_WATERPLANE_TO_VOLUME_MAX = 12.67

DSYHS_S3_PRISMATIC_MIN = DSYHS_S4_PRISMATIC_MIN = DSYHS_2008_PRISMATIC_MIN = 0.52
DSYHS_S3_PRISMATIC_MAX = DSYHS_S4_PRISMATIC_MAX = DSYHS_2008_PRISMATIC_MAX = 0.60

DSYHS_S3_LWL_TO_VOLUME_MIN = DSYHS_S4_LWL_TO_VOLUME_MIN = DSYHS_2008_LWL_TO_VOLUME_MIN = 4.34
DSYHS_S3_LWL_TO_VOLUME_MAX = DSYHS_S4_LWL_TO_VOLUME_MAX = DSYHS_2008_LWL_TO_VOLUME_MAX = 8.5

DSYHS_S4_LWL_TO_BWL_MIN = DSYHS_2008_LWL_TO_BWL_MIN = 2.73
DSYHS_S3_LWL_TO_BWL_MIN = 2.76
DSYHS_2008_LWL_TO_BWL_MAX = 5.88
DSYHS_S3_LWL_TO_BWL_MAX = DSYHS_S4_LWL_TO_BWL_MAX = 5.0

DSYHS_S3_LCB_FRACTION_MIN = DSYHS_S4_LCB_FRACTION_MIN = DSYHS_2008_LCB_FRACTION_MIN = 0.5
DSYHS_2008_LCB_FRACTION_MAX = DSYHS_S4_LCB_FRACTION_MAX = 0.582
DSYHS_S3_LCB_FRACTION_MAX = 0.56

DSYHS_S4_LCF_FRACTION_MIN = DSYHS_2008_LCF_FRACTION_MIN = 0.518
DSYHS_S4_LCF_FRACTION_MAX = DSYHS_2008_LCF_FRACTION_MAX = 0.595

DSYHS_2008_MIDSHIP_MIN = 0.65
DSYHS_2008_MIDSHIP_MAX = 0.79

check_bwl_to_draft_2008 = partial(check_domain, domain=DSYHS_2008, quantity="bwl/Tc", lower=DSYHS_2008_BEAM_TO_DRAFT_MIN, upper=DSYHS_2008_BEAM_TO_DRAFT_MAX)
check_bwl_to_draft_s4 = partial(check_domain, domain=DSYHS_S4, quantity="bwl/Tc", lower=DSYHS_S4_BEAM_TO_DRAFT_MIN, upper=DSYHS_S4_BEAM_TO_DRAFT_MAX)
check_bwl_to_draft_s3 = partial(check_domain, domain=DSYHS_S3, quantity="bwl/Tc", lower=DSYHS_S3_BEAM_TO_DRAFT_MIN, upper=DSYHS_S3_BEAM_TO_DRAFT_MAX)
check_bwl_to_draft = check_bwl_to_draft_s3

check_waterplane_to_volume_s4 = partial(check_domain, domain=DSYHS_S4, quantity="Aw / Vc**(2./3.)", lower=DSYHS_S4_WATERPLANE_TO_VOLUME_MIN, upper=DSYHS_S4_WATERPLANE_TO_VOLUME_MAX)

check_prismatic_2008 = partial(check_domain, domain=DSYHS_2008, quantity="Cp", lower=DSYHS_2008_PRISMATIC_MIN, upper=DSYHS_2008_PRISMATIC_MAX)
check_prismatic_s4 = partial(check_domain, domain=DSYHS_S4, quantity="Cp", lower=DSYHS_S4_PRISMATIC_MIN, upper=DSYHS_S4_PRISMATIC_MAX)
check_prismatic_s3 = partial(check_domain, domain=DSYHS_S3, quantity="Cp", lower=DSYHS_S3_PRISMATIC_MIN, upper=DSYHS_S3_PRISMATIC_MAX)

check_lwl_to_volume_2008 = partial(check_domain, domain=DSYHS_2008, quantity="lwl / Vc**(1./3.)", lower=DSYHS_2008_LWL_TO_VOLUME_MIN, upper=DSYHS_2008_LWL_TO_VOLUME_MAX)
check_lwl_to_volume_s4 = partial(check_domain, domain=DSYHS_S4, quantity="lwl / Vc**(1./3.)", lower=DSYHS_S4_LWL_TO_VOLUME_MIN, upper=DSYHS_S4_LWL_TO_VOLUME_MAX)
check_lwl_to_volume_s3 = partial(check_domain, domain=DSYHS_S3, quantity="lwl / Vc**(1./3.)", lower=DSYHS_S3_LWL_TO_VOLUME_MIN, upper=DSYHS_S3_LWL_TO_VOLUME_MAX)
check_lwl_to_volume = check_lwl_to_volume_s3

check_lwl_to_bwl_2008 = partial(check_domain, domain=DSYHS_2008, quantity="lwl / bwl", lower=DSYHS_2008_LWL_TO_BWL_MIN, upper=DSYHS_2008_LWL_TO_BWL_MAX)
check_lwl_to_bwl_s4 = partial(check_domain, domain=DSYHS_S4, quantity="lwl / bwl", lower=DSYHS_S4_LWL_TO_BWL_MIN, upper=DSYHS_S4_LWL_TO_BWL_MAX)
check_lwl_to_bwl_s3 = partial(check_domain, domain=DSYHS_S3, quantity="lwl / bwl", lower=DSYHS_S3_LWL_TO_BWL_MIN, upper=DSYHS_S3_LWL_TO_BWL_MAX)
check_lwl_to_bwl = check_lwl_to_bwl_s3

check_lcb_fraction_2008 = partial(check_domain, domain=DSYHS_2008, quantity="LCB fraction", lower=DSYHS_2008_LCB_FRACTION_MIN, upper=DSYHS_2008_LCB_FRACTION_MAX)
check_lcb_fraction_s4 = partial(check_domain, domain=DSYHS_S4, quantity="LCB fraction", lower=DSYHS_S4_LCB_FRACTION_MIN, upper=DSYHS_S4_LCB_FRACTION_MAX)
check_lcb_fraction_s3 = partial(check_domain, domain=DSYHS_S3, quantity="LCB fraction", lower=DSYHS_S3_LCB_FRACTION_MIN, upper=DSYHS_S3_LCB_FRACTION_MAX)
check_lcb_fraction = check_lcb_fraction_s3

check_lcf_fraction_2008 = partial(check_domain, domain=DSYHS_2008, quantity="LCF fraction", lower=DSYHS_2008_LCF_FRACTION_MIN, upper=DSYHS_2008_LCF_FRACTION_MAX)
check_lcf_fraction_s4 = partial(check_domain, domain=DSYHS_S4, quantity="LCF fraction", lower=DSYHS_S4_LCF_FRACTION_MIN, upper=DSYHS_S4_LCF_FRACTION_MAX)

check_midship_2008 = partial(check_domain, domain=DSYHS_2008, quantity="Cm", lower=DSYHS_2008_MIDSHIP_MIN, upper=DSYHS_2008_MIDSHIP_MAX)


# Appendages resistance functions


def keel_residuary_delta_heeled_ks(boatspeed: float,
                                   heel_angle: float,
                                   Vc: float,
                                   Vk: float,
                                   Tc: float,
                                   T: float,
                                   lwl: float,
                                   bwl: float,
                                   coe: Union[Tuple[float, float, float], None] = None,
                                   rho_water: float = RHO_SEA_WATER_20C) -> Force:
    """ Calculates the keel residuary resistance [N] delta at a given heel angle

    KS stands for Keuning Sonnenberg.
    Keuning and Sonnenberg wrote "Approximation of the Hydrodynamic Forces
    on a Sailing Yacht based on the 'Delft Systematic Yacht Hull Series'" in 1998.
    keel_residuary_resistance_delta_heeled_KS uses the results from this publication.

    boatspeed : boat speed [m/s]
    heel_angle : The heel_angle [degrees], between -90 and 90.
    Vc : Canoe body volume of displaced water [m**3] , must be > 0.
    Vk : Keel volume of displaced water [m**3] , must be > 0.
    Tc : Canoe body draft [m], must be >= 0
    T : Total draft [m], must be >= 0
    lwl : Length at Waterline [m] , must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    coe : tuple of 3 floats, optional (default [lwl / 2, 0, -Tc/3]).
        Represents the x, y, z coordinates of the centre of effort.
        For residuary resistance, the surfacic barycentre of the underwater
        hull is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    rho_water :  Water density [kg/m**3], must be >= 0 (default RHO_SEA_WATER)

    Returns a Force object, representing the residuary resistance [N] delta of the keel at any heel angle.

    """
    if Vc <= 0:
        raise ValueError("Vc should be strictly positive")
    if Vk <= 0:
        raise ValueError("Vk should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if T <= 0:
        raise ValueError("T should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    # DSYHS checks
    # midship area coefficient cannot be checked here
    check_lwl_to_bwl(value=lwl / bwl)
    check_bwl_to_draft(value=bwl / Tc)
    check_lwl_to_volume(value=lwl / Vc**(1./3.))

    if coe is None:
        coe = (lwl / 2., 0., -Tc / 3.)

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.

    froude = froude_number(speed=boatspeed, lwl=lwl)

    rrk_delta_at_heel_angle = (Vk * rho_water * GRAVITY_STANDARD) * \
                              (-3.5837 * Tc / T +
                               -0.0518 * bwl / Tc +
                               0.5958 * bwl / T +
                               0.2055 * lwl / Vc ** (1. / 3.)) \
                              * froude ** 2 * radians(abs(heel_angle))

    return Force(-rrk_delta_at_heel_angle * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])


@memoize
def _get_keel_residuary_ks_interpolator(Vc: float,
                                        Vk: float,
                                        Tc: float,
                                        T: float,
                                        bwl: float,
                                        Zcbk: float,
                                        rho_water: float) -> interp1d:
    froude_numbers_table = [0.00, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60]
    A0 = [0.00000, -0.00104, -0.00550, -0.01110, -0.00713, -0.03581, -0.00470, 0.00553, 0.04822, 0.01021]
    A1 = [0.00000, 0.00172, 0.00597, 0.01421, 0.02632, 0.08649, 0.11592, 0.07371, 0.00660, 0.14173]
    A2 = [0.00000, 0.00117, 0.00390, 0.00069, -0.00232, 0.00999, -0.00064, 0.05991, 0.07048, 0.06409]
    A3 = [0.00000, -0.00008, -0.00009, 0.00021, 0.00039, 0.00017, 0.00035, -0.00114, -0.00035, -0.00192]

    # (A0[i] + A1[i]*T/Bwl + A2[i]*Math.pow(Tc+Zcbk,3)/Vk + A3[i]*Vc/Vk);
    rrk = [(Vk * rho_water * GRAVITY_STANDARD) * (A0[i] + A1[i] * T / bwl
                                                  + A2[i] * (Tc + Zcbk) / Vk ** (1.00 / 3.00)
                                                  + A3[i] * Vc / Vk)
           for i in range(len(froude_numbers_table))]
    return interp1d(froude_numbers_table, rrk, kind='slinear', bounds_error=True)


def keel_residuary_ks(boatspeed: float,
                      Vc: float,
                      Vk: float,
                      Tc: float,
                      T: float,
                      lwl: float,
                      bwl: float,
                      Zcbk: float,
                      coe: Union[Tuple[float, float, float], None] = None,
                      rho_water: float = RHO_SEA_WATER_20C) -> Force:
    """Calculate the residuary resistance [N] of the keel at 0° heel.

    KS stands for Keuning Sonnenberg.
    Keuning and Sonnenberg wrote "Approximation of the Hydrodynamic Forces
    on a Sailing Yacht based on the 'Delft Systematic Yacht Hull Series'"
    in 1998.
    keel_residuary_resistance_KS uses the results from this publication.

    boatspeed : boat speed [m/s]
    Vc : Canoe body volume of displaced water [m**3] , must be > 0.
    Vk : Keel volume of displaced water [m**3] , must be > 0.
    Tc : Canoe body draft [m], must be >= 0
    T : Total draft [m], must be >= 0
    lwl : Length at Waterline [m] , must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    Zcbk : Distance to the center of buoyancy of the keel [m]n from top of keel. Must be > 0.
    coe : tuple of 3 floats, optional (default (lwl / 2, 0, -Tc/3)).
        Represents the x, y, z coordinates of the centre of effort.
        For residuary resistance, the surfacic barycentre of the underwater hull is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    rho_water : Water density [kg/m**3], must be >= 0 (default RHO_SEA_WATER)

    Returns a Force object, representing the residuary resistance [N] of the keel at 0° heel.

    """
    if Vc <= 0:
        raise ValueError("Vc should be strictly positive")
    if Vk <= 0:
        raise ValueError("Vk should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if T <= 0:
        raise ValueError("T should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")
    # DSYHS checks
    # midship area coefficient cannot be checked here
    check_lwl_to_bwl(value=lwl / bwl)
    check_bwl_to_draft(value=bwl / Tc)
    check_lwl_to_volume(value=lwl / Vc**(1./3.))

    if coe is None:
        coe = (lwl / 2., 0., -Tc / 3.) if coe is None else coe

    interp = _get_keel_residuary_ks_interpolator(Vc, Vk, Tc, T, bwl, Zcbk, rho_water)

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.

    froude = froude_number(speed=boatspeed, lwl=lwl)
    if froude <= 0.6:
        residuary_resistance = max([0.,
                                   interp(abs(boatspeed / sqrt(GRAVITY_STANDARD * lwl)))])
    else:
        at_0_59 = interp(0.59)
        at_0_60 = interp(0.60)
        residuary_resistance = at_0_60 + \
                               (froude - 0.6) * (at_0_60 - at_0_59) / 0.01

    return Force(-residuary_resistance * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])


# Hull resistance functions

# **** CHECKLIST ****
# Tests / Checks against known values where possible              OK
#     100.79 instead of 99 ??


def due_to_heel(boatspeed: float,
                heel_angle: float,
                lwl: float,
                bwl: float,
                Sc: float,
                Tc: float,
                T: float,
                coe: Union[Tuple[float, float, float], None] = None,
                rho_water: float = RHO_SEA_WATER_20C) -> Force:
    """ Calculate the drag [N] due to heel.
    LA stands for Larsson, one of the authors of Principles of Yacht Design.
    The calculation is described on page 83 of the book Principles of Yacht Design

    boatspeed : boat speed [m/s]
    heel_angle : The heel_angle [degrees], between -90 and 90.
    lwl : Length at Waterline [m] , must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    Sc : WSAc - Canoe body Wetted Surface Area [m**2], must be >= 0
    Tc : Canoe body draft [m], must be >= 0
    T : Total draft [m], must be >= 0
    coe : tuple of 3 floats, optional (default [lwl / 2, 0, -Tc/3]).
        Represents the x, y, z coordinates of the centre of effort.
        For residuary resistance, the surfacic barycentre of the underwater hull is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    rho_water : Water density [kg/m**3], must be >= 0

    Returns a Force object, representing the drag [N] due to heel.

    References
    ----------
    http://klz-resistance.blogspot.fr/

    """
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Sc <= 0:
        raise ValueError("Sc should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if T <= 0:
        raise ValueError("T should be strictly positive")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    if coe is None:
        coe = (lwl / 2., 0., -Tc / 3.)

    Ch = ((6.747 * (Tc / T)) +
          (2.517 * (bwl / Tc)) +
          (3.710 * (bwl / Tc) * (Tc / T))) * 0.001

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.
    froude = abs(boatspeed / sqrt(GRAVITY_STANDARD * lwl))

    heel_resistance = (0.5 * rho_water * boatspeed ** 2 * Sc * Ch * froude ** 2 * radians(abs(heel_angle)))

    return Force(-heel_resistance * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])


@memoize
def _get_hull_residuary_delta_heeled_ks_interpolator(Vc: float,
                                                     lwl: float,
                                                     bwl: float,
                                                     Tc: float,
                                                     lcb_fraction: float,
                                                     rho_water: float) -> interp1d:

    froude_numbers_table = [0.00, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55]

    u0 = [0.0000, -0.0268, 0.6628, 1.6433, -0.8659, -3.2715, -0.1976, 1.5873]
    u1 = [0.0000, -0.0014, -0.0632, -0.2144, -0.0354, 0.1372, -0.1480, -0.3749]
    u2 = [0.0000, -0.0057, -0.0699, -0.1640, 0.2226, 0.5547, -0.6593, -0.7105]
    u3 = [0.0000, 0.0016, 0.0069, 0.0199, 0.0188, 0.0268, 0.1862, 0.2146]
    u4 = [0.0000, -0.0070, 0.0459, -0.0540, -0.5800, -1.0064, -0.7489, -0.4818]
    u5 = [0.0000, -0.0017, -0.0004, -0.0268, -0.1133, -0.2026, -0.1648, -0.1174]

    # Conversion du LCB format 0.5xy en LCB format -x.y
    LCB = (0.5 - lcb_fraction) * 100

    delta_rh_at_20_degree_heel = [(Vc * rho_water * GRAVITY_STANDARD) *
                                  (0.001 * u0[i] + 0.001 * u1[i] * lwl / bwl + 0.001 * u2[
                                      i] * bwl / Tc
                                   + 0.001 * u3[i] * (bwl / Tc) ** 2 + 0.001 * u4[i] * LCB + 0.001 * u5[
                                       i] * LCB ** 2)
                                  for i in range(len(froude_numbers_table))]

    return interp1d(froude_numbers_table, delta_rh_at_20_degree_heel, kind='slinear', bounds_error=True)


def hull_residuary_delta_heeled_ks(boatspeed: float,
                                   heel_angle: float,
                                   Vc: float,
                                   lwl: float,
                                   bwl: float,
                                   Tc: float,
                                   lcb_fraction: float = 0.53,
                                   coe: Union[Tuple[float, float, float], None] = None,
                                   rho_water: float = RHO_SEA_WATER_20C) -> Force:
    """Calculate the hull residuary resistance [N] delta at a given heel angle.
    KS stands for Keuning Sonnenberg.
    Keuning and Sonnenberg wrote "Approximation of the Hydrodynamic Forces
    on a Sailing Yacht based on the 'Delft Systematic Yacht Hull Series'"
    in 1998.
    hull_residuary_resistance_delta_heeled_KS uses the results from this publication.

    boatspeed : boat speed [m/s]
    heel_angle : The heel_angle [degrees], between -90 and 90.
    Vc : Canoe body volume of displaced water [m**3] , must be > 0.
    lwl : Length at Waterline [m] , must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    Tc : Canoe body draft [m], must be >= 0
    lcb_fraction : Longitudinal position of Center of Buoyancy from bow.
        It is in fraction format: e.g. 0.560. Must be > 0 and < 1.
        The DSYHS range is 0.500 to 0.582.
    coe : tuple of 3 floats, optional (default [lwl / 2, 0, -Tc / 3])
        Represents the x, y, z coordinates of the centre of effort.
        For residuary resistance, the surfacic barycentre of the underwater hull is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    rho_water : float, optional (default RHO_SEA_WATER)
        Water density [kg/m**3], must be >= 0.

    """
    if Vc <= 0:
        raise ValueError("Vc should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if not 0. < lcb_fraction < 1.:
        raise ValueError("lcb_fraction should be between 0 and 1")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    # DSYHS checks
    # midship area coefficient cannot be checked here
    check_lwl_to_bwl(value=lwl / bwl)
    check_bwl_to_draft(value=bwl / Tc)
    check_lwl_to_volume(value=lwl / Vc**(1./3.))
    check_lcb_fraction(value=lcb_fraction)

    coe = (lwl / 2., 0., -Tc / 3.) if coe is None else coe

    interp = _get_hull_residuary_delta_heeled_ks_interpolator(Vc, lwl, bwl, Tc, lcb_fraction, rho_water)

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.

    froude = froude_number(speed=boatspeed, lwl=lwl)
    if froude <= 0.55:
        rrh_delta_at_20_deg_heel = max([0., interp(froude)])
        rrh_delta_at_heel_angle = (rrh_delta_at_20_deg_heel *
                                   6.0 *
                                   (radians(abs(heel_angle))) ** 1.7)
    else:
        at_0_54 = interp(0.54)
        at_0_55 = interp(0.55)
        rrh_delta_at_20_deg_heel = (at_0_55 +
                                    (froude - 0.55) *
                                    (at_0_55 - at_0_54) / 0.01)
        rrh_delta_at_heel_angle = (rrh_delta_at_20_deg_heel * 6.0 *
                                   (radians(abs(heel_angle))) ** 1.7)

    return Force(-rrh_delta_at_heel_angle * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])


@memoize
def _get_hull_residuary_resistance_2008_interpolator(Vc: float,
                                                     lwl: float,
                                                     bwl: float,
                                                     Tc : float,
                                                     Aw: float,
                                                     Cm: float,
                                                     lcb_fraction: float,
                                                     lcf_fraction: float,
                                                     Cp: float,
                                                     rho_water: float) -> interp1d:
    froude_numbers_table = [0.00, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75]

    a0 = [0.0000, -0.0005, -0.0003, -0.0002, -0.0009, -0.0026, -0.0064, -0.0218, -0.0388, -0.0347, -0.0361,  0.0008,  0.0108,  0.1023]
    a1 = [0.0000,  0.0023,  0.0059, -0.0156,  0.0016, -0.0567, -0.4034, -0.5261, -0.5986, -0.4764,  0.0037,  0.3728, -0.1238,  0.7726]
    a2 = [0.0000, -0.0086, -0.0064,  0.0031,  0.0337,  0.0446, -0.1250, -0.2945, -0.3038, -0.2361, -0.2960, -0.3667, -0.2026,  0.5040]
    a3 = [0.0000, -0.0015,  0.0070, -0.0021, -0.0285, -0.1091,  0.0273,  0.2485,  0.6033,  0.8726,  0.9661,  1.3957,  1.1282,  1.7867]
    a4 = [0.0000,  0.0061,  0.0014, -0.0070, -0.0367, -0.0707, -0.1341, -0.2428, -0.0430,  0.4219,  0.6123,  1.0343,  1.1836,  2.1934]
    a5 = [0.0000,  0.0010,  0.0013,  0.0148,  0.0218,  0.0914,  0.3578,  0.6293,  0.8332,  0.8990,  0.7534,  0.3230,  0.4973, -1.5479]
    a6 = [0.0000,  0.0001,  0.0005,  0.0010,  0.0015,  0.0021,  0.0045,  0.0081,  0.0106,  0.0096,  0.0100,  0.0072,  0.0038, -0.0115]
    a7 = [0.0000,  0.0052, -0.0020, -0.0043, -0.0172, -0.0078,  0.1115,  0.2086,  0.1336, -0.2272, -0.3352, -0.4632, -0.4477, -0.0977]

    Rrh = [Vc * rho_water * GRAVITY_STANDARD *
           (a0[i]
            + (a1[i] * lcb_fraction + a2[i] * Cp + a3[i] * (Vc ** 0.6666666666) / Aw + a4[i] * bwl / lwl) * ((Vc ** 0.3333333333) / lwl)
            + (a5[i] * lcb_fraction / lcf_fraction + a6[i] * bwl / Tc + a7[i] * Cm) * ((Vc ** 0.3333333333) / lwl))
           for i in range(len(froude_numbers_table))]

    # TODO : slinear ??
    return interp1d(froude_numbers_table, Rrh, kind='slinear', bounds_error=True)


def hull_residuary_resistance_2008(boatspeed: float,
                                   Vc: float,
                                   lwl: float,
                                   bwl: float,
                                   Tc: float,
                                   Aw: float,
                                   Cm: float,
                                   lcb_fraction: float = 0.53,
                                   lcf_fraction: float = 0.56,
                                   Cp: float = 0.56,
                                   coe: Union[Tuple[float, float, float], None] = None,
                                   rho_water: float = RHO_SEA_WATER_20C) -> Force:
    """Calculate the residuary resistance [N] on the hull
    (source: p24 of Fossati's Aero-hydrodynamics ...)

    boatspeed : boat speed [m/s]
    Vc : Canoe body volume of displaced water [m**3], must be > 0.
    lwl : Length at Waterline [m], must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    Aw : Waterplane area [m**2] , must be > 0.
    Cm : Midship section coefficient, must be >= 0
    lcb_fraction : Longitudinal position of Center of Buoyancy from bow.
        It is in fraction format e.g. 0.560. Must be > 0 and < 1. The DSYHS range is 0.500 to 0.582.
    lcf_fraction : Longitudinal position of Center of Flotation from bow.
        It is in fraction format e.g. 0.590. Must be > 0 and < 1. The DSYHS range is 0.518 to 0.595.
    Cp : Prismatic Coefficient, Must be > 0 and < 1. The DSYHS range is 0.520 to 0.600.
    coe : tuple of 3 floats, optional (default [lwl / 2, 0, 0])
        Represents the x, y, z coordinates of the centre of effort.
        For residuary resistance, the surfacic barycentre of the underwater hull is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    rho_water : Water density [kg/m**3], must be >= 0.

    References
    ----------
    Aero-hydrodynamics and the performance of sailing yachts. p24

    """
    if Vc <= 0:
        raise ValueError("Vc should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Aw <= 0:
        raise ValueError("Aw should be strictly positive")
    if Cm <= 0:
        raise ValueError("Cm should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if not 0. < lcb_fraction < 1.:
        raise ValueError("lcb_fraction should be between 0 and 1")
    if not 0. < lcf_fraction < 1.:
        raise ValueError("lcf_fraction should be between 0 and 1")
    if not 0. < Cp < 1.:
        raise ValueError("Cp should be between 0 and 1")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    # DSYHS checks
    # bwl / Tc cannot be checked here
    # midship area coefficient cannot be checked here
    check_lwl_to_bwl_2008(value=lwl / bwl)
    check_lwl_to_volume_2008(value=lwl / Vc**(1./3.))
    check_lcb_fraction_2008(value=lcb_fraction)
    check_lcf_fraction_2008(value=lcf_fraction)
    check_prismatic_2008(value=Cp)
    check_midship_2008(value=Cm)

    coe = (lwl / 2., 0., 0.) if coe is None else coe

    interp = _get_hull_residuary_resistance_2008_interpolator(
        Vc, lwl, bwl, Tc, Aw, Cm, lcb_fraction, lcf_fraction, Cp, rho_water)

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.

    froude = froude_number(speed=boatspeed, lwl=lwl)

    if froude <= 0.75:
        residuary_resistance = max([0., interp(froude)])
    else:
        # Rough extrapolation using the derivative at the end of the domain
        at_0_74 = interp(0.74)
        at_0_75 = interp(0.75)
        residuary_resistance = at_0_75 + (froude - 0.75) * (at_0_75 - at_0_74) / 0.01

    return Force(-residuary_resistance * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])


# **** CHECKLIST ****
# Tests / Checks against know values where possible              OK !!
#     525 instead of 611


@memoize
def _get_hull_residuary_resistance_ks_series4_interpolator(Vc: float,
                                                           lwl: float,
                                                           bwl: float,
                                                           Aw: float,
                                                           Sc: float,
                                                           lcb_fraction: float,
                                                           lcf_fraction: float,
                                                           Cp: float,
                                                           rho_water: float) -> interp1d:
    froude_numbers_table = [0.00, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60]

    a0 = [0.0000, -0.0014, 0.0004, 0.0014, 0.0027, 0.0056, 0.0032, -0.0064, -0.0171, -0.0201, 0.0495, 0.0808]
    a1 = [0.0000, 0.0403, -0.1808, -0.1071, 0.0463, -0.8005, -0.1011, 2.3095, 3.4017, 7.1576, 1.5618, -5.3233]
    a2 = [0.0000, 0.0470, 0.1793, 0.0637, -0.1263, 0.4891, -0.0813, -1.5152, -1.9862, -6.3304, -6.0661, -1.1513]
    a3 = [0.0000, -0.0227, -0.0004, 0.0090, 0.0150, 0.0269, -0.0382, 0.0751, 0.3242, 0.5829, 0.8641, 0.9663]
    a4 = [0.0000, -0.0119, 0.0097, 0.0153, 0.0274, 0.0519, 0.0320, -0.0858, -0.1450, 0.1630, 1.1702, 1.6084]
    a5 = [0.0000, 0.0061, 0.0118, 0.0011, -0.0299, -0.0313, -0.1481, -0.5349, -0.8043, -0.3966, 1.7610, 2.7459]
    a6 = [0.0000, -0.0086, -0.0055, 0.0012, 0.0110, 0.0292, 0.0837, 0.1715, 0.2952, 0.5023, 0.9176, 0.8491]
    a7 = [0.0000, -0.0307, 0.1721, 0.1021, -0.0595, 0.7314, 0.0223, -2.4550, -3.5284, -7.1579, -2.1191, 4.7129]
    a8 = [0.0000, -0.0553, -0.1728, -0.0648, 0.1220, -0.3619, 0.1587, 1.1865, 1.3575, 5.2534, 5.4281, 1.1089]

    Rrh = [Vc * rho_water * GRAVITY_STANDARD *
           (a0[i] + (a1[i] * lcb_fraction + a2[i] * Cp + a3[i] * (Vc ** 0.6666666666) / Aw
                     + a4[i] * bwl / lwl) * ((Vc ** 0.3333333333) / lwl)
            + (a5[i] * (Vc ** 0.6666666666) / Sc + a6[i] * lcb_fraction / lcf_fraction
               + a7[i] * (lcb_fraction ** 2) + a8[i] * (Cp ** 2)) * ((Vc ** 0.3333333333) / lwl))
           for i in range(len(froude_numbers_table))]

    return interp1d(froude_numbers_table, Rrh, kind='slinear', bounds_error=True)


def hull_residuary_resistance_ks_series4(boatspeed: float,
                                         Vc: float,
                                         lwl: float,
                                         bwl: float,
                                         Aw: float,
                                         Sc: float,
                                         lcb_fraction: float = 0.53,
                                         lcf_fraction: float = 0.56,
                                         Cp: float = 0.56,
                                         coe: Union[Tuple[float, float, float], None] = None,
                                         rho_water: float = RHO_SEA_WATER_20C) -> Force:
    """Calculate the residuary resistance [N] on the hull, according to the DSYHS Series 4 regression
    KS stands for Keuning Sonnenberg.
    Keuning and Sonnenberg wrote "Approximation of the Hydrodynamic Forces
    on a Sailing Yacht based on the 'Delft Systematic Yacht Hull Series'"
    in 1998. hull_residuary_resistance_KS uses the results from this publication

    boatspeed : boat speed [m/s]
    Vc : Canoe body volume of displaced water [m**3], must be > 0.
    lwl : Length at Waterline [m], must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    Aw : Waterplane area [m**2] , must be > 0.
    Sc : WSAc - Canoe body Wetted Surface Area [m**2], must be >= 0
    lcb_fraction : Longitudinal position of Center of Buoyancy from bow.
        It is in fraction format e.g. 0.560. Must be > 0 and < 1. The DSYHS range is 0.500 to 0.582.
    lcf_fraction : Longitudinal position of Center of Flotation from bow.
        It is in fraction format e.g. 0.590. Must be > 0 and < 1. The DSYHS range is 0.518 to 0.595.
    Cp : Prismatic Coefficient, Must be > 0 and < 1. The DSYHS range is 0.520 to 0.600.
    coe : tuple of 3 floats, optional (default [lwl / 2, 0, 0])
        Represents the x, y, z coordinates of the centre of effort.
        For residuary resistance, the surfacic barycentre of the underwater hull is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    rho_water : Water density [kg/m**3], must be >= 0.

    References
    ----------
    Approximation of the Hydrodynamic Forces on a Sailing Yacht
    based on the 'Delft Systematic Yacht Hull Series (by Keuning & Sonnenberg)

    http://radiosailingtechnology.com/index.php/hulls/
                           drag-measurements-on-an-international-one-metre-yacht

    http://klz-resistance.blogspot.fr/

    """
    if Vc <= 0:
        raise ValueError("Vc should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Aw <= 0:
        raise ValueError("Aw should be strictly positive")
    if Sc <= 0:
        raise ValueError("Sc should be strictly positive")
    if not 0. < lcb_fraction < 1.:
        raise ValueError("lcb_fraction should be between 0 and 1")
    if not 0. < lcf_fraction < 1.:
        raise ValueError("lcf_fraction should be between 0 and 1")
    if not 0. < Cp < 1.:
        raise ValueError("Cp should be between 0 and 1")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    # DSYHS checks
    # bwl / Tc cannot be checked here
    # midship area coefficient cannot be checked here
    check_lwl_to_bwl_s4(value=lwl / bwl)
    check_lwl_to_volume_s4(value=lwl / Vc**(1./3.))
    check_lcb_fraction_s4(value=lcb_fraction)
    check_lcf_fraction_s4(value=lcf_fraction)
    check_prismatic_s4(value=Cp)
    check_waterplane_to_volume_s4(value=Aw / Vc**(2./3.))

    coe = (lwl / 2., 0., 0.) if coe is None else coe

    interp = _get_hull_residuary_resistance_ks_series4_interpolator(
        Vc, lwl, bwl, Aw,  Sc, lcb_fraction, lcf_fraction, Cp, rho_water)

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.

    froude = froude_number(speed=boatspeed, lwl=lwl)
    if froude <= 0.6:
        residuary_resistance = max([0., interp(froude)])
    else:
        # Rough extrapolation using the derivative at the end of the domain
        at_0_59 = interp(0.59)
        at_0_60 = interp(0.60)
        residuary_resistance = at_0_60 + (froude - 0.6) * (at_0_60 - at_0_59) / 0.01

    return Force(-residuary_resistance * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])


@memoize
def _get_hull_residuary_resistance_la_series3_interpolator(Vc: float,
                                                           lwl: float,
                                                           bwl: float,
                                                           Tc: float,
                                                           Aw: float,
                                                           lcb_fraction: float,
                                                           Cp: float,
                                                           rho_water: float) -> interp1d:
    lcb_percent = (0.5 - lcb_fraction) * 100.

    froude_numbers_table = [0.00, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250,
                            0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425,
                            0.450, 0.475, 0.500, 0.525, 0.550, 0.575, 0.600,
                            0.625, 0.650, 0.675, 0.700, 0.725, 0.750]

    a0 = [0.000000, -6.735654, -0.382870, -1.503526, 11.29218, 22.17867, 25.90867, 40.97559, 45.83759, 89.20382,
          212.6788, 336.2354, 566.5476, 743.4107, 1200.620, 180.1004, 243.9994, 282.9873, 313.4109, 337.0038,
          356.4572, 324.7357, 301.1268, 292.0571, 284.4641, 256.6367, 304.1803]
    a1 = [0.000000, 38.36831, 38.17290, 24.40803, -14.51947, -49.16784, -74.75668, -114.2855, -184.7646, -393.0127,
          -801.7908, -1085.134, -1609.632, -1708.263, -2751.715, -31.50257, -44.52551, -51.51953, -56.58257,
          -59.19029, -62.85395, -51.31252, -39.79631, -31.85303, -25.14558, -19.31922, -30.11512]
    a2 = [0.000000, -0.008193, 0.007243, 0.012200, 0.047182, 0.085998, 0.153521, 0.207226, 0.357031, 0.617466,
          1.087307, 1.644191, 2.016090, 2.435809, 3.208577, -7.451141, -11.15456, -12.97310, -14.41978, -16.06975,
          -16.85112, -15.34595, -15.02299, -15.58548, -16.15423, -13.08450, -15.85429]
    a3 = [0.000000, 0.055234, 0.026644, 0.067221, 0.085176, 0.150725, 0.188568, 0.250827, 0.338343, 0.460472,
          0.538938, 0.532702, 0.265722, 0.013553, 0.254920, 2.195042, 2.179046, 2.274505, 2.326117, 2.419156,
          2.437056, 2.334146, 2.059657, 1.847926, 1.703981, 2.152824, 2.863173]
    a4 = [0.000000, -1.997242, -5.295332, -2.448582, -2.673016, -2.878684, -0.889467, -3.072662, 3.871658, 11.54327,
          10.80273, -1.224173, -29.24412, -81.16189, -132.0424, 2.689623, 3.857403, 4.343662, 4.690432, 4.766793,
          5.078768, 3.855368, 2.545676, 1.569917, 0.817912, 0.348305, 1.524379]
    a5 = [0.000000, -38.86081, -39.55032, -31.91370, -11.41819, 7.167049, 24.12137, 53.01570, 132.2568, 331.1197,
          667.6445, 831.1445, 1154.091, 937.4014, 1489.269, 0.006480, 0.009676, 0.011066, 0.012147, 0.014147,
          0.014980, 0.013695, 0.013588, 0.014014, 0.014575, 0.011343, 0.014031]
    a6 = [0.000000, 0.956591, 1.219563, 2.216098, 5.654065, 8.600272, 10.48516, 13.02177, 10.86054, 8.598136,
          12.39815, 26.18321, 51.46175, 115.6006, 196.3406]
    a7 = [0.000000, -0.002171, 0.000052, 0.000074, 0.007021, 0.012981, 0.025348, 0.035934, 0.066809, 0.104073,
          0.166473, 0.238795, 0.288046, 0.365071, 0.528225]
    a8 = [0.000000, 0.272895, 0.824568, 0.244345, -0.094934, -0.327085, -0.854940, -0.715457, -1.719215, -2.815203,
          -3.026131, -2.450470, -0.178354, 1.838967, 1.379102]
    a9 = [0.000000, -0.017516, -0.047842, -0.015887, 0.006325, 0.018271, 0.048449, 0.039874, 0.095977, 0.155960,
          0.165055, 0.139154, 0.018446, -0.062023, 0.013577]

    Rrh = np.zeros(27)

    for i, froude in enumerate(froude_numbers_table):
        if froude <= 0.451:
            Rrh[i] = Vc * rho_water * GRAVITY_STANDARD * 0.001 * \
                     (
                         a0[i]
                         + a1[i] * Cp
                         + a2[i] * lcb_percent
                         + a3[i] * (bwl / Tc)
                         + a4[i] * (lwl / Vc ** (1.0 / 3.0))
                         + a5[i] * Cp ** 2.0
                         + a6[i] * (Cp * lwl / Vc ** (1.0 / 3.0))
                         + a7[i] * lcb_percent ** 2.0
                         + a8[i] * (lwl / Vc ** (1.0 / 3.0)) ** 2.0
                         + a9[i] * (lwl / Vc ** (1.0 / 3.0)) ** 3.0)
        else:
            Rrh[i] = Vc * rho_water * GRAVITY_STANDARD * 0.001 * \
                     (a0[i]
                      + a1[i] * (lwl / bwl)
                      + a2[i] * (Aw / Vc ** (2.0 / 3.0))
                      + a3[i] * lcb_percent
                      + a4[i] * (lwl / bwl) ** 2.0
                      + a5[i] * (lwl / bwl) * ((Aw / Vc ** (2.0 / 3.0)) ** 3.0))

    return interp1d(froude_numbers_table, Rrh, kind='slinear', bounds_error=True)


def hull_residuary_resistance_la_series3(boatspeed: float,
                                         Vc: float,
                                         lwl: float,
                                         bwl: float,
                                         Tc: float,
                                         Aw: float,
                                         lcb_fraction: float = 0.53,
                                         Cp: float = 0.56,
                                         coe: Union[Tuple[float, float, float], None] = None,
                                         rho_water: float = RHO_SEA_WATER_20C) -> Force:
    """The hull upright residuary resistance (N) according to the DSYHS Series 3 regression

    Calculates the residuary resistance [N] on the hull.
    LA stands for Larsson, one of the authors of Principles of Yacht Design.
    The calculation is described on pages 75 and 76 of the book
    Principles of Yacht Design

    boatspeed : boat speed [m/s]
    Vc : Canoe body volume of displaced water [m**3] , must be > 0.
    lwl : Length at Waterline [m] , must be > 0.
    bwl : Beam at Waterline [m] , must be > 0.
    Tc : Canoe body draft [m], must be >= 0
    Aw : Waterplane area [m**2] , must be > 0.
    lcb_fraction : Longitudinal position of Center of Buoyancy from bow.
        It is in fraction format: e.g. 0.560. Must be > 0 and < 1.
        The DSYHS range is 0.500 to 0.582.
    Cp : Prismatic Coefficient. Must be > 0 and < 1. The DSYHS range is 0.520 to 0.600.
    coe : tuple of 3 floats, optional (default [lwl / 2, 0, -Tc / 3]).
        Represents the x, y, z coordinates of the centre of effort.
        For residuary resistance, the surfacic barycentre of the underwater
        hull is an acceptable estimate.
        x=lwl / 2 , y = 0., z = -Tc/3 is another acceptable estimate.
    rho_water : Water density [kg/m**3], must be >= 0.

    Returns a Force object, representing the residuary resistance [N] on the hull.

    References
    ----------
    Principles of Yacht Design - Larsson & Eliasson
    http://radiosailingtechnology.com/index.php/hulls/
                           drag-measurements-on-an-international-one-metre-yacht
    http://klz-resistance.blogspot.fr/

    """
    if Vc <= 0:
        raise ValueError("Vc should be strictly positive")
    if lwl <= 0:
        raise ValueError("lwl should be strictly positive")
    if bwl <= 0:
        raise ValueError("bwl should be strictly positive")
    if Tc <= 0:
        raise ValueError("Tc should be strictly positive")
    if Aw <= 0:
        raise ValueError("Aw should be strictly positive")
    if not 0. < lcb_fraction < 1.:
        raise ValueError("lcb_fraction should be between 0 and 1")
    if not 0. < Cp < 1.:
        raise ValueError("Cp should be between 0 and 1")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    # DSYHS checks (cf page 74 of principles of yacht design
    check_lwl_to_bwl_s3(value=lwl / bwl)
    check_bwl_to_draft_s3(value=bwl / Tc)
    check_lwl_to_volume_s3(value=lwl / Vc**(1./3.))
    check_lcb_fraction_s3(value=lcb_fraction)
    check_prismatic_s3(value=Cp)

    coe = (lwl / 2., 0., -Tc / 3.) if coe is None else coe

    interp = _get_hull_residuary_resistance_la_series3_interpolator(
        Vc, lwl, bwl, Tc, Aw, lcb_fraction, Cp, rho_water)

    boatspeed_sign = boatspeed / abs(boatspeed) if boatspeed != 0. else 0.
    froude = froude_number(speed=boatspeed, lwl=lwl)
    if froude <= 0.75:
        residuary_resistance = max([0., interp(froude)])
    else:
        at_0_74 = interp(0.74)
        at_0_75 = interp(0.75)
        residuary_resistance = (at_0_75 +
                                (froude - 0.75) *
                                (at_0_75 - at_0_74) / 0.01)

    return Force(-residuary_resistance * boatspeed_sign, 0., 0., coe[0], coe[1], coe[2])
