# coding: utf-8

r"""Savitsky method for planning hulls"""

from typing import Callable, Union, Tuple
import logging
from math import sqrt, radians, cos, tan, sin
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from matplotlib import pyplot
from ydeos_hydrodynamics.constants import RHO_SEA_WATER_20C, GRAVITY_STANDARD, \
    KINEMATIC_VISCOSITY_WATER_20C
from ydeos_hydrodynamics.friction_lines import ittc57
from ydeos_hydrodynamics.water import RHO_WATER_MIN, RHO_WATER_MAX
from ydeos_hydrodynamics.force import Force
from ydeos_hydrodynamics.memoize import memoize

logger = logging.getLogger(__name__)


def root_scalar_no_fallback(func: Callable,
                            x0: float,
                            x1: float,
                            xtol: float,
                            min_value: float) -> Union[float, None]:
    sol = root_scalar(func, method="secant", x0=x0, x1=x1, xtol=xtol)
    if sol.converged is True and type(sol.root) == float and sol.root > min_value:
        return sol.root
    else:
        return None


def root_scalar_with_fallback(func: Callable,
                              x0: float,
                              x1: float,
                              bracket: Tuple[float, float],
                              xtol: float,
                              min_value: float) -> Union[float, None]:
    r"""Root finding that tries another method if the first did not work properly"""
    sol = root_scalar(func, method="secant", x0=x0, x1=x1, xtol=xtol)
    if sol.converged is True and type(sol.root) == float and sol.root > min_value:
        return sol.root
    else:
        sol = root_scalar(func, method="brentq", bracket=bracket, xtol=xtol)
        if sol.converged is True and type(sol.root) == float and sol.root > min_value:
            return sol.root
        else:
            return None


def _cl_beta_cl_zero(boatspeed: float,
                     m: float,
                     b: float,
                     beta: float,
                     rho_water: float = RHO_SEA_WATER_20C,
                     gravity: float = GRAVITY_STANDARD) -> Tuple[float, float]:
    r"""Determine the lift coefficients for a 'beta' deadrise angle (-> Cl_beta)
    and for a zero deadrise angle (-> Cl_zero)

    boatspeed : [m/s]
    m : Mass displacement [kg]
    b : maximum beam between chines (or between spray rails) [m]
    beta : deadrise angle (take the average of the angles at the transom and at the CG) [degrees]
           zero for an absolutely flat bottom

    Returns a tuple of (Cl_beta, Cl_zero)

    """
    # Lift coefficient for beta deadrise angle
    Cl_beta = m * gravity / (0.5 * rho_water * boatspeed**2 * b**2)

    # The check might seem useless but floating point precision near zero may play bad tricks !!
    if type(Cl_beta) == float and Cl_beta > 0:
        # Determine Cl_zero
        def func(cl_z):
            r"""Function for which we search a root"""
            # Cl_beta = cl_zero - 0.0065 * beta * cl_zero**0.6
            return cl_z - 0.0065 * beta * cl_z**0.6 - Cl_beta

        Cl_zero = root_scalar_with_fallback(func, x0=Cl_beta, x1=10*Cl_beta,
                                            bracket=(0., 100.), xtol=1e-6, min_value=0.)

        if Cl_zero is not None:
            return Cl_beta, Cl_zero
        else:
            raise RuntimeError("Could not find a solution for Cl_zero")
    else:
        raise RuntimeError("Cl_beta is expected to be a strictly positive float, but isn't")


def _lambda(trim_angle: float, Cv: float, Cl_zero: float) -> float:
    r"""Determine lambda (mean wetted length to beam ratio)

    trim_angle [degrees] must be positive (exponentiation by 1.1 would create
    a complex number if trim_angle < 0)
    Cv : Beam Froude Number / Speed coefficient
    Cl_zero : lift coefficient for a zero deadrise angle
    """
    #
    def func(lambd):
        r"""Function for which we search a root"""
        return trim_angle**1.1 * (0.012 * lambd**0.5 + 0.0055 * lambd**2.5 / Cv**2) - Cl_zero

    # Could not use the fallback version as the bracket elements have the same sign
    lambda_ = root_scalar_no_fallback(func, x0=0., x1=10., xtol=1e-6, min_value=0.)

    if lambda_ is not None:
        return lambda_
    else:
        # It is ugly, yet better than finding no solution
        lambdas = np.linspace(0, 100, 1000)
        results = [func(l) for l in lambdas]
        if min(results) > 0. or max(results) < 0.:
            # no chance to find a solution
            # TODO : is it better to raise or to return an arbitrary realistic value ?
            # raise RuntimeError("Could not find a solution for lambda")
            return 1.
        else:
            i = interp1d(results, lambdas)
            return i(0.)


def bow_down_moment(boatspeed: float,
                    trim_angle: float,
                    m: float,
                    LCG: float,
                    VCG: float,
                    b: float,
                    epsilon: float,
                    beta: float,
                    f: float,
                    rho_water: float = RHO_SEA_WATER_20C,
                    kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
                    gravity: float = GRAVITY_STANDARD) -> Tuple[float, ...]:
    r"""Computes the bow down moment

    boatspeed : [m/s]
    trim_angle : bow up is positive [degrees]
    m : Mass displacement [kg]
    LCG : Distance from transom to CG [m]
    VCG : distance from baseline (keel) to CG [m]
    b : maximum beam between chines (or between spray rails) [m]
    epsilon : propeller shaft inclination relative to baseline [degrees]
    beta : deadrise angle (take the average of the angles at the transom and at the CG) [degrees]
           zero for an absolutely flat bottom
    f : distance between shaftline and CG [m]

    Returns the bow-down moment [N.m]

    """
    Cv = boatspeed / sqrt(gravity * b)  # Speed coefficient / Beam Froude Number
    Cl_beta, Cl_zero = _cl_beta_cl_zero(boatspeed, m, b, beta, rho_water, gravity)
    lambda_ = _lambda(trim_angle, Cv, Cl_zero)
    Lm = lambda_ * b  # Mean wetted length
    Cf = ittc57(boatspeed, Lm, kinematic_viscosity)

    # Delta lambda (Increase in length-beam ratio due to spray)
    # TODO : implement the values from the graph on page 189 of Principles of Yacht Design
    delta_lambda = 0.

    # Frictional resistance
    beta_rad = radians(beta)
    Rf = Cf * 0.5 * rho_water * boatspeed**2 * (lambda_ + delta_lambda) * b**2 / cos(beta_rad)

    ff = VCG - b / 4. * tan(beta_rad)  # ff : Lever arm for Rf [m]

    # Skipped steps 11 and 12 (appendages Ra and fa
    # TODO : implement the appendages in the savitsky method
    Ra = 0
    fa = 0

    # Lcp : Distance from centre of pressure to trailing edge
    Lcp = Lm * (0.75 - (1 / (5.21 * Cv**2 / lambda_**2 + 2.39)))

    # e : lever arm for the pressure force
    e = LCG - Lcp

    # M : bow down moment
    epsilon_rad = radians(epsilon)
    trim_angle_rad = radians(trim_angle)
    Mh = gravity * m * (e * cos(trim_angle_rad + epsilon_rad) / cos(epsilon_rad) - f * sin(trim_angle_rad) / cos(epsilon_rad))
    Mf = Rf * (ff - e * tan(epsilon_rad) - f / cos(epsilon_rad))
    Ma = Ra * (fa - e * tan(epsilon_rad) - f / cos(epsilon_rad))
    M = Mh + Mf + Ma

    R = (gravity * m * sin(trim_angle_rad) + Rf) * cos(trim_angle_rad + epsilon_rad) / cos(epsilon_rad)

    return M, Mh, Mf, Ma, R, Rf, Ra, ff, fa, Cv, Cl_beta, Cl_zero, lambda_, delta_lambda, Lcp, e


@memoize
def savitsky(boatspeed: float,
             m: float,
             LCG: float,
             VCG: float,
             b: float,
             epsilon: float,
             beta: float,
             f: float,
             rho_water: float = RHO_SEA_WATER_20C,
             kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
             gravity: float = GRAVITY_STANDARD) -> Tuple[float, ...]:
    r"""Savitsky method for planning hull.

    WARNING : Appendages and delta lambda NOT IMPLEMENTED

    Example values on p204 example in Principles of Yacht Design

    The Savitsky method is iterative and needs to find an equilibrium
    trim angle (i.e. bow down moment = 0) before computing the resistance.
    The procedure is explained in Principles of Yacht Design p185-195 (procedure p192)

    trim angles must remains positive as the exponentiation by 1.1 in _lambda would create a complex

    boatspeed : [m/s]
    m : Mass displacement [kg]
    LCG : Distance from transom to CG [m]
    VCG : distance from baseline (keel) to CG [m]
    b : maximum beam between chines (or between spray rails) [m]
    epsilon : propeller shaft inclination relative to baseline [degrees]
    beta : deadrise angle (take the average of the angles at the transom and at the CG) [degrees]
           zero for an absolutely flat bottom
    f : distance between shaftline and CG [m]

    """
    if boatspeed <= 0:
        raise ValueError("boatspeed should be strictly positive")
    if m <= 0:
        raise ValueError("m (mass) should be strictly positive")
    if LCG <= 0:
        raise ValueError("LCG should be strictly positive")
    if VCG <= 0:
        raise ValueError("VCG should be strictly positive")
    if b <= 0:
        raise ValueError("b (beam between chines) should be strictly positive")
    if not -45. <= epsilon <= 45.:
        raise ValueError("epsilon (propeller shaft angle) is unrealistic")
    if not 0. <= beta <= 45.:
        raise ValueError("beta (deadrise angle) is unrealistic")
    if not RHO_WATER_MIN < rho_water < RHO_WATER_MAX:
        raise ValueError(f"rho_water should be between {RHO_WATER_MIN} and {RHO_WATER_MAX}")

    def func(trim_angle):
        r"""The function for which we search a root"""
        return bow_down_moment(boatspeed, trim_angle, m, LCG, VCG, b, epsilon,
                               beta, f, rho_water, kinematic_viscosity, gravity)[0]

    tau_0M = root_scalar_with_fallback(func, x0=0., x1=45.,
                                       bracket=(0., 45), xtol=1e-6, min_value=0.)

    if tau_0M is not None:
        zero_bow_down_moment_trim_angle = tau_0M
        logger.debug("Solver")
    else:
        trim_angles = np.linspace(0., 45., 1000)
        Ms = [bow_down_moment(boatspeed,
                              trim_angle,
                              m,
                              LCG,
                              VCG,
                              b,
                              epsilon,
                              beta,
                              f,
                              rho_water,
                              kinematic_viscosity,
                              gravity)[0] for trim_angle in trim_angles]

        i = interp1d(list(Ms), trim_angles)

        zero_bow_down_moment_trim_angle = i(0.)
        logger.debug("Brute force interpolation")

    # Recompute everything using the equilibrium trim angle
    M, Mh, Mf, Ma, R, Rf, Ra, ff, fa, Cv, Cl_beta, Cl_zero, lambda_, delta_lambda, Lcp, e = \
        bow_down_moment(boatspeed,
                        zero_bow_down_moment_trim_angle,
                        m,
                        LCG,
                        VCG,
                        b,
                        epsilon,
                        beta,
                        f,
                        rho_water,
                        kinematic_viscosity,
                        gravity)

    return (float(zero_bow_down_moment_trim_angle),
            M, Mh, Mf, Ma,
            R, Rf, Ra,
            ff, fa,
            Cv, Cl_beta, Cl_zero, lambda_, delta_lambda, Lcp, e)


def savitsky_force(boatspeed: float,
                   m: float,
                   LCG: float,
                   VCG: float,
                   b: float,
                   epsilon: float,
                   beta: float,
                   f: float,
                   rho_water: float = RHO_SEA_WATER_20C,
                   kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
                   gravity: float = GRAVITY_STANDARD) -> Force:
    r"""Savitsky output formatted as a Force"""
    _, _, _, _, _, R, _, _, _, _, _, _, _, _, _, Lcp, _ = \
        savitsky(boatspeed, m, LCG, VCG, b, epsilon, beta, f, rho_water, kinematic_viscosity, gravity)
    return Force(-R, 0., 0., Lcp, 0, 0)


def savitsky_plot_moment_vs_trim_angle(boatspeed: float,
                                       m: float,
                                       LCG: float,
                                       VCG: float,
                                       b: float,
                                       epsilon: float,
                                       beta: float,
                                       f: float,
                                       rho_water: float = RHO_SEA_WATER_20C,
                                       kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
                                       gravity: float = GRAVITY_STANDARD) -> None:
    r"""Plot the moments vs trim angles.
    This procedure is mostly aimed at debugging in case of problem"""
    trim_angles = np.linspace(0., 45., 1000)
    Ms = [bow_down_moment(boatspeed,
                          trim_angle,
                          m,
                          LCG,
                          VCG,
                          b,
                          epsilon,
                          beta,
                          f,
                          rho_water,
                          kinematic_viscosity,
                          gravity)[0] for trim_angle in trim_angles]

    pyplot.scatter(trim_angles, Ms, marker='+')
    pyplot.show()
    pyplot.scatter(Ms, trim_angles, marker='+')
    pyplot.show()


def savitsky_report(boatspeed: float,
                    m: float,
                    LCG: float,
                    VCG: float,
                    b: float,
                    epsilon: float,
                    beta: float,
                    f: float,
                    rho_water: float = RHO_SEA_WATER_20C,
                    kinematic_viscosity: float = KINEMATIC_VISCOSITY_WATER_20C,
                    gravity: float = GRAVITY_STANDARD) -> None:
    r"""Print a report to stdout"""
    zero_bow_down_moment_trim_angle, M, Mh, Mf, Ma, R, Rf, Ra, ff, fa, Cv, Cl_beta, Cl_zero, lambda_, delta_lambda, Lcp, e = \
        savitsky(boatspeed, m, LCG, VCG, b, epsilon, beta, f, rho_water, kinematic_viscosity, gravity)
    print("  Trim angle : %.3f" % zero_bow_down_moment_trim_angle)
    print("           M : %.3f" % M)
    print("          Mh : %.3f" % Mh)
    print("          Mf : %.3f" % Mf)
    print("          Ma : %.3f" % Ma)
    print("           R : %.3f" % R)
    print("          Rf : %.3f" % Rf)
    print("          Ra : %.3f" % Ra)
    print("          ff : %.3f" % ff)
    print("          fa : %.3f" % fa)
    print("          Cv : %.3f" % Cv)
    print("     Cl_beta : %.3f" % Cl_beta)
    print("     Cl_zero : %.3f" % Cl_zero)
    print("      lambda : %.3f" % lambda_)
    print("delta_lambda : %.3f" % delta_lambda)
    print("         Lcp : %.3f" % Lcp)
    print("           e : %.3f" % e)
