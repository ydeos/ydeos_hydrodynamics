#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Benchmarking of the hydrodynamics functions"""

from ydeos_benchmark.benchmark import run_benchmark_simple

from ydeos_hydrodynamics.holtrop_mennen import holtrop_mennen
from ydeos_hydrodynamics.savitski import savitsky
from ydeos_hydrodynamics.geometry import shift_coe_around_x
from ydeos_hydrodynamics.water import water_fresh_density, water_sea_density, \
    water_fresh_dynamic_viscosity, water_sea_dynamic_viscosity, \
    water_fresh_kinematic_viscosity, water_sea_kinematic_viscosity
from ydeos_hydrodynamics.delft import hull_residuary_resistance_2008, \
    due_to_heel, hull_residuary_delta_heeled_ks, hull_residuary_resistance_la_series3, \
    hull_residuary_resistance_ks_series4, keel_residuary_delta_heeled_ks, \
    keel_residuary_ks
from ydeos_hydrodynamics.mit import hull_residuary_resistance_mit
from ydeos_hydrodynamics.friction_lines import appendage_friction_coefficient, \
    bulb_friction_coefficient, hughes, ittc57, laminar, prandtl, transition_x
from ydeos_hydrodynamics.form_factors import hoerner_streamlined_body, \
    hoerner, k_watanabe
from ydeos_hydrodynamics.viscous import appendage_viscous, ballast_viscous, \
    hull_viscous
from ydeos_hydrodynamics.sideforce import beta_zero, keel_downwash_on_rudder, \
    lift_curve_slope, sideforce_by_element_fossati, sideforce_production_keuning_verwerft
from ydeos_hydrodynamics.orc import create_interpolation, hull_residuary_resistance_orc_2013
from ydeos_hydrodynamics.froude import froude_number, speed_ms
from ydeos_hydrodynamics.reynolds import reynolds_number
from ydeos_hydrodynamics.righting_moment import rm_estimate_gerritsma, rm_orc2013, \
    trm_estimate
from ydeos_hydrodynamics.wetted_area import wetted_area


# List of functions to benchmark
to_profile = (
    (holtrop_mennen, [25.*(1852/3600), 37500, 205., 32., 10., 0.5833, 0.5075, "U", 7381.45, 16., 0.98, 0.45, 0.75,
                      20., 4.0, 10., (50.,), (200.,), (.5,)], {}),
    (savitsky, [10, 1000, 3, 0.5, 3, 10, 10, 1.], {}),
    (shift_coe_around_x, [(1., 2., 3.), 90.], {}),
    (water_fresh_density, [20.], {}),
    (water_sea_density, [20.], {}),
    (water_fresh_dynamic_viscosity, [20.], {}),
    (water_sea_dynamic_viscosity, [20.], {}),
    (water_fresh_kinematic_viscosity, [20.], {}),
    (water_sea_kinematic_viscosity, [20.], {}),
    (hull_residuary_resistance_2008, [], {'boatspeed': 1., 'Vc': 0.04, 'lwl': 1., 'bwl': 0.2, 'Tc': 0.06, 'Aw': 0.1,
                                          'Cm': 0.75, 'lcb_fraction': 0.52, 'lcf_fraction': 0.53, 'Cp': 0.55}),
    #
    (appendage_friction_coefficient, [], {'speed': 1., 'dimension': 1., 'thickness_to_chord': 0.06}),
    (appendage_viscous, [], {'boatspeed': 1., 'heel_angle': 10., 'coe': (1., 2., 3.), 'average_chord': 0.1,
                             'thickness_to_chord': 0.06, 'wetted_area': 0.01}),
    (ballast_viscous, [], {'boatspeed': 1., 'heel_angle': 10., 'coe': (1., 2., 3.), 'bulb_length': 0.3,
                           'bulb_diameter': 0.05, 'bulb_wetted_area': 0.02}),
    (beta_zero, [], {'beam_to_draft_ratio': 4., 'heel_angle_degrees': 20.}),
    (bulb_friction_coefficient, [], {'speed': 1., 'dimension': 0.4}),
    # (check_btr, [3.], {}),
    # (check_bwl_to_draft, [3., 1.], {}),
    # (check_lcb_fraction, [0.5], {}),
    # (check_lcf_fraction, [0.5], {}),
    # (check_lcf_fraction_s4, [0.5], {}),
    # (check_lvr, [10.], {}),
    # (check_lwl_to_bwl, [5., 1.], {}),
    # (check_lwl_to_volume, [10., 2.], {}),
    # (check_lwl_to_volume_s3, [10., 2.], {}),
    # (check_lwl_to_volume_s4, [10., 2.], {}),
    # (check_prismatic, [0.55], {}),
    # (check_prismatic_s3, [0.55], {}),
    # (check_prismatic_s4, [0.55], {}),
    # (check_waterplane_to_volume, [42., 10.], {}),
    # (check_waterplane_to_volume_s4, [42., 10.], {}),
    (create_interpolation, [], {'Xs': [x for x in range(40)], 'Ys': [y for y in range(40)],
                                'Zs': [z for z in range(40*40)]}),
    (due_to_heel, [], {'boatspeed': 1., 'heel_angle': 10., 'lwl': 1., 'bwl': 0.2, 'Sc': 0.1, 'Tc': 0.06, 'T': 0.4, 'coe': None}),
    (froude_number, [1., 1.], {}),
    (hoerner, [0.05], {}),
    (hoerner_streamlined_body, [0.05], {}),
    (hughes, [], {'speed': 1., 'dimension': 1.}),
    (hull_residuary_delta_heeled_ks, [], {'boatspeed': 1., 'heel_angle': 10., 'Vc': 0.04, 'lwl': 1., 'bwl': 0.2,
                                          'Tc': 0.06, 'lcb_fraction': 0.53}),
    (hull_residuary_resistance_ks_series4, [], {'boatspeed': 1., 'Vc': 0.04, 'lwl': 1., 'bwl': 0.2, 'Aw': 0.1,
                                                'Sc': 0.12, 'lcb_fraction': 0.52, 'lcf_fraction':0.53, 'Cp': 0.55}),
    (hull_residuary_resistance_la_series3, [], {'boatspeed': 1., 'Vc': 0.04, 'lwl': 1., 'bwl': 0.2, 'Tc': 0.06,
                                                'Aw': 0.1, 'lcb_fraction': 0.52, 'Cp': 0.55}),
    (hull_residuary_resistance_mit, [], {'boatspeed': 1., 'Vc': 0.04, 'lwl': 1., 'bwl': 0.2, 'Tc': 0.06}),
    (hull_residuary_resistance_orc_2013, [], {'boatspeed': 1., 'Vc': 0.004, 'lwl': 1., 'bwl': 0.2, 'Tc': 0.055}),
    (hull_viscous, [], {'boatspeed': 1., 'lwl': 1., 'Sc': 0.12, 'coe': (1., 2., 3.)}),
    # (is_within_limits_btr, [], {'BTR': 22.}),
    # (is_within_limits_bwl_to_draft, [], {'bwl': 3., 'Tc': 1., 'lower': 2., 'upper': 4.}),
    # (is_within_limits_lcb_fraction, [], {'lcb_fraction':0.5, 'lower': 0.4, 'upper': 0.6}),
    # (is_within_limits_lcf_fraction, [], {'lcf_fraction': 0.5}),
    # (is_within_limits_lcf_fraction_s4, [], {'lcf_fraction': 0.5}),
    # (is_within_limits_lvr, [], {'LVR': 42.}),
    # (is_within_limits_lwl_to_bwl, [], {'lwl': 1., 'bwl': 0.2, 'lower': 4., 'upper': 6.}),
    # (is_within_limits_lwl_to_volume, [], {'lwl': 1., 'Vc': 0.04}),
    # (is_within_limits_lwl_to_volume_s3, [], {'lwl': 1., 'Vc': 0.04}),
    # (is_within_limits_lwl_to_volume_s4, [], {'lwl': 1., 'Vc': 0.04}),
    # (is_within_limits_prismatic, [], {'Cp': 0.55}),
    # (is_within_limits_prismatic_s3, [], {'Cp': 0.55}),
    # (is_within_limits_prismatic_s4, [], {'Cp': 0.55}),
    # (is_within_limits_waterplane_to_volume, [], {'Aw': 0.1, 'Vc': 0.04}),
    # (is_within_limits_waterplane_to_volume_s4, [], {'Aw': 0.1, 'Vc': 0.04}),
    (ittc57, [], {'speed': 1., 'dimension': 1.}),
    (k_watanabe, [], {'lwl': 1., 'bwl': 0.2, 'Tc': 0.06, 'block_coefficient': 0.4}),
    (keel_downwash_on_rudder, [], {'keel_Cl': 0.4, 'keel_effective_aspect_ratio': 3., 'heel_angle_degrees': 20.}),
    (keel_residuary_delta_heeled_ks, [], {'boatspeed': 1., 'heel_angle': 10., 'Vc': 0.04, 'Vk': 0.001, 'Tc': 0.06,
                                          'T': 0.4, 'lwl': 1., 'bwl': 0.2, 'coe': (1., 2., 3.)}),
    (keel_residuary_ks, [], {'boatspeed': 1., 'Vc': 0.04, 'Vk': 0.01, 'Tc': 0.06, 'T': 0.4, 'lwl': 1., 'bwl': 0.2,
                             'Zcbk': 0.2, 'coe': (1., 2., 3.)}),
    (laminar, [], {'speed': 1., 'dimension': 1.}),
    (lift_curve_slope, [], {'effective_aspect_ratio': 3., 'sweep_back_angle': 10.}),
    (prandtl, [], {'speed': 1., 'dimension': 1.}),
    (reynolds_number, [], {'speed': 1., 'dimension': 1.}),
    (rm_estimate_gerritsma, [], {'boatspeed': 1., 'heel_angle': 10., 'displacement': 4., 'lwl': 1., 'bwl': 0.2,
                                 'Tc': 0.06, 'Gdwl': 0.2}),
    (rm_orc2013, [], {'heel_angle': 10., 'displacement': 4., 'lwl': 1., 'bwl': 0.2, 'Tc': 0.06}),
    (sideforce_by_element_fossati, [], {'boatspeed': 1., 'heel_angle': 10., 'leeway_angle': 4., 'rudder_angle': 1.,
                                        'coe': (1., 2., 3.), 'projected_area': 0.05, 'span': 0.4}),
    (sideforce_production_keuning_verwerft, [], {'boatspeed': 1., 'heel_angle': 10., 'leeway_angle': 4.,
                                                 'rudder_angle': 1., 'keel_projected_area': 0.02, 'keel_span': 0.4,
                                                 'rudder_projected_area': 0.01, 'rudder_span': 0.2, 'lwl': 1.,
                                                 'bwl': 0.2, 'Tc': 0.06, 'hull_coe': (1., 2., 3.),
                                                 'keel_coe': (1., 2., 3.), 'rudder_coe': (1., 2., 3.),
                                                 'keel_sweep_back_angle': 10., 'rudder_sweep_back_angle': 0.}),
    (speed_ms, [], {'froude': 0.4, 'lwl': 1.}),
    (transition_x, [], {'speed': 1., 'transition_reynolds_number': 1000000.}),
    (trm_estimate, [], {'trim_angle': 2., 'displacement': 4., 'lwl': 1., 'bwl': 0.2, 'Tc': 0.06, 'Gdwl': 0.2}),
    (wetted_area, [], {'upright_wsa': 0.1, 'bwl': 0.2, 'Tc': 0.06, 'Cm': 0.73, 'heel': 20.}),)


if __name__ == "__main__":
    n = 1000
    run_benchmark_simple(to_profile, n_times=n)
    # run_benchmark_complete(to_profile, n_times=n, save_results=True)
