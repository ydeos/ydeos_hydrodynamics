#!/usr/bin/env python
# coding: utf-8

r"""Tests for the sideforce.py module."""

import pytest
from ydeos_hydrodynamics.sideforce import sideforce_by_element_fossati, \
    sideforce_production_keuning_verwerft
from ydeos_forces.forces import Force, SystemOfForces


spbef = sideforce_by_element_fossati(boatspeed=4.,
                                     heel_angle=0.,
                                     leeway_angle=4.,
                                     rudder_angle=0.,
                                     coe=(0.5, 0., -0.1),
                                     projected_area=0.03,
                                     span=0.2,
                                     aspect_ratio_multiplier=2.,
                                     use_cos_heel_angle=True,
                                     rudder_angle_influence=0.5,
                                     rho_water=1025.)


def test_known_value_sideforce_by_element_fossati():
    r"""No know data."""
    assert True


def test_zero_boatspeed_sideforce_by_element_fossati():
    r"""Zero boatspeed."""
    f = sideforce_by_element_fossati(boatspeed=0.,
                                     heel_angle=0.,
                                     leeway_angle=4.,
                                     rudder_angle=0.,
                                     coe=(0.5, 0., -0.1),
                                     projected_area=0.03,
                                     span=0.2,
                                     aspect_ratio_multiplier=2.,
                                     use_cos_heel_angle=True,
                                     rudder_angle_influence=0.5,
                                     rho_water=1025.)
    assert f.fx == 0.
    assert f.fy == 0.
    assert f.fz == 0.


def test_negative_projected_area_sideforce_by_element_fossati():
    r"""Negative projected area."""
    with pytest.raises(ValueError):
        sideforce_by_element_fossati(boatspeed=4.,
                                     heel_angle=0.,
                                     leeway_angle=4.,
                                     rudder_angle=0.,
                                     coe=(0.5, 0., -0.1),
                                     projected_area=-0.03,
                                     span=0.2,
                                     aspect_ratio_multiplier=2.,
                                     use_cos_heel_angle=True,
                                     rudder_angle_influence=0.5,
                                     rho_water=1025.)


def test_sign_and_symmetry_sideforce_by_element_fossati():
    r"""Test signs and symmetry."""
    force = sideforce_by_element_fossati(boatspeed=4.,
                                         heel_angle=0.,
                                         leeway_angle=4.,
                                         rudder_angle=0.,
                                         coe=(0.5, 0., -0.1),
                                         projected_area=0.03,
                                         span=0.2,
                                         aspect_ratio_multiplier=2.,
                                         use_cos_heel_angle=True,
                                         rudder_angle_influence=0.5,
                                         rho_water=1025.)
    force_backwards = sideforce_by_element_fossati(boatspeed=-4.,
                                                   heel_angle=0.,
                                                   leeway_angle=4.,
                                                   rudder_angle=0.,
                                                   coe=(0.5, 0., -0.1),
                                                   projected_area=0.03,
                                                   span=0.2,
                                                   aspect_ratio_multiplier=2.,
                                                   use_cos_heel_angle=True,
                                                   rudder_angle_influence=0.5,
                                                   rho_water=1025.)
    force_nose_down = sideforce_by_element_fossati(boatspeed=4.,
                                                   heel_angle=0.,
                                                   leeway_angle=-4.,
                                                   rudder_angle=0.,
                                                   coe=(0.5, 0., -0.1),
                                                   projected_area=0.03,
                                                   span=0.2,
                                                   aspect_ratio_multiplier=2.,
                                                   use_cos_heel_angle=True,
                                                   rudder_angle_influence=0.5,
                                                   rho_water=1025.)
    assert force.fx == -force_backwards.fx
    assert force.fy == -force_nose_down.fy
    assert force.fx < 0


def test_heel_angle_resistance_sideforce_by_element_fossati():
    r"""x component of force with opposite heel angles should be equal."""
    force = sideforce_by_element_fossati(boatspeed=4.,
                                         heel_angle=20.,
                                         leeway_angle=4.,
                                         rudder_angle=3.,
                                         coe=(0.5, 0., -0.1),
                                         projected_area=0.03,
                                         span=0.2,
                                         aspect_ratio_multiplier=2.,
                                         use_cos_heel_angle=True,
                                         rudder_angle_influence=0.5,
                                         rho_water=1025.)
    force_heeled_to_windward_negative_leeway = \
        sideforce_by_element_fossati(boatspeed=4.,
                                     heel_angle=-20.,
                                     leeway_angle=-4.,
                                     rudder_angle=-3.,
                                     coe=(0.5, 0., -0.1),
                                     projected_area=0.03,
                                     span=0.2,
                                     aspect_ratio_multiplier=2.,
                                     use_cos_heel_angle=True,
                                     rudder_angle_influence=0.5,
                                     rho_water=1025.)
    assert force.fx == force_heeled_to_windward_negative_leeway.fx


def test_heel_angle_sideforce_sideforce_by_element_fossati():
    r"""Symmetric situations, mirrored by XZ plane."""
    force = sideforce_by_element_fossati(boatspeed=4.,
                                         heel_angle=20.,
                                         leeway_angle=4.,
                                         rudder_angle=3.,
                                         coe=(0.5, 0., -0.1),
                                         projected_area=0.03,
                                         span=0.2,
                                         aspect_ratio_multiplier=2.,
                                         use_cos_heel_angle=True,
                                         rudder_angle_influence=0.5,
                                         rho_water=1025.)
    force_heeled_to_windward_negative_leeway = \
        sideforce_by_element_fossati(boatspeed=4.,
                                     heel_angle=-20.,
                                     leeway_angle=-4.,
                                     rudder_angle=-3.,
                                     coe=(0.5, 0., -0.1),
                                     projected_area=0.03,
                                     span=0.2,
                                     aspect_ratio_multiplier=2.,
                                     use_cos_heel_angle=True,
                                     rudder_angle_influence=0.5,
                                     rho_water=1025.)
    assert force.fy == -force_heeled_to_windward_negative_leeway.fy


def test_rudder_angle_sideforce_by_element_fossati():
    r"""Opposite rudder angles."""
    force = sideforce_by_element_fossati(boatspeed=4.,
                                         heel_angle=0.,
                                         leeway_angle=0.,
                                         rudder_angle=5.,
                                         coe=(0.5, 0., -0.1),
                                         projected_area=0.03,
                                         span=0.2,
                                         aspect_ratio_multiplier=2.,
                                         use_cos_heel_angle=True,
                                         rudder_angle_influence=0.5,
                                         rho_water=1025.)
    force_negative_rudder = \
        sideforce_by_element_fossati(boatspeed=4.,
                                     heel_angle=0.,
                                     leeway_angle=0.,
                                     rudder_angle=-5.,
                                     coe=(0.5, 0., -0.1),
                                     projected_area=0.03,
                                     span=0.2,
                                     aspect_ratio_multiplier=2.,
                                     use_cos_heel_angle=True,
                                     rudder_angle_influence=0.5,
                                     rho_water=1025.)
    assert force.fx == force_negative_rudder.fx
    assert force.fy == -force_negative_rudder.fy


def test_coe_position_sideforce_by_element_fossati():
    r"""Centre of effort position tests."""
    force = sideforce_by_element_fossati(boatspeed=4.,
                                         heel_angle=20.,
                                         leeway_angle=5.,
                                         rudder_angle=5.,
                                         coe=(0.5, 0., -0.1),
                                         projected_area=0.03,
                                         span=0.2,
                                         aspect_ratio_multiplier=2.,
                                         use_cos_heel_angle=True,
                                         rudder_angle_influence=0.5,
                                         rho_water=1025.)
    assert force.px == 0.5
    assert force.py < 0.
    assert force.pz > - 0.2167

    force = sideforce_by_element_fossati(boatspeed=4.,
                                         heel_angle=-20.,
                                         leeway_angle=5.,
                                         rudder_angle=5.,
                                         coe=(0.5, 0., -0.1),
                                         projected_area=0.03,
                                         span=0.2,
                                         aspect_ratio_multiplier=2.,
                                         use_cos_heel_angle=True,
                                         rudder_angle_influence=0.5,
                                         rho_water=1025.)
    assert force.px == 0.5
    assert force.py > 0.
    assert force.pz > - 0.2167


def test_speed_multiplier_sideforce_by_element_fossati():
    r"""Different speed multipliers."""
    force = sideforce_by_element_fossati(boatspeed=4.,
                                         heel_angle=20.,
                                         leeway_angle=5.,
                                         rudder_angle=5.,
                                         coe=(0.5, 0., -0.1),
                                         projected_area=0.03,
                                         span=0.2,
                                         aspect_ratio_multiplier=2.,
                                         use_cos_heel_angle=True,
                                         rudder_angle_influence=0.5,
                                         speed_multiplier=1.0)
    # force = spbef.force(boatspeed=4.,
    #                     heel_angle=20., leeway_angle=5., rudder_angle=5.)
    force_if_mult = sideforce_by_element_fossati(boatspeed=4.,
                                                 heel_angle=20.,
                                                 leeway_angle=5.,
                                                 rudder_angle=5.,
                                                 coe=(0.5, 0., -0.1),
                                                 projected_area=0.03,
                                                 span=0.2,
                                                 aspect_ratio_multiplier=2.,
                                                 use_cos_heel_angle=True,
                                                 rudder_angle_influence=0.5,
                                                 speed_multiplier=1.2)
    assert force_if_mult.fx / 1.2 ** 2 - 1e-6 <= force.fx <= force_if_mult.fx / 1.2 ** 2 + 1e-6
    assert force_if_mult.fy / 1.2 ** 2 - 1e-6 <= force.fy <= force_if_mult.fy / 1.2 ** 2 + 1e-6


lwl = 1.0
bwl = 0.15893
Tc = 0.05525
bulb_length = 0.4
bulb_diameter = bulb_length * 0.089169
spkv09 = \
    sideforce_production_keuning_verwerft(boatspeed=4.,
                                          heel_angle=0.,
                                          leeway_angle=4.,
                                          rudder_angle=0.,
                                          keel_projected_area=0.025240,
                                          keel_span=0.42 - Tc - bulb_diameter,
                                          rudder_projected_area=0.012402,
                                          rudder_span=0.248,
                                          lwl=lwl,
                                          bwl=bwl,
                                          Tc=Tc,
                                          # hull_coe=None,
                                          keel_coe=(0.460, 0., - 0.2167),
                                          rudder_coe=(0.0536, 0., -0.11495))


spkv09_swept_20 = \
    sideforce_production_keuning_verwerft(boatspeed=4.,
                                          heel_angle=0.,
                                          leeway_angle=4.,
                                          rudder_angle=0.,
                                          keel_projected_area=0.025240,
                                          keel_span=0.42 - Tc - bulb_diameter,
                                          rudder_projected_area=0.012402,
                                          rudder_span=0.248,
                                          lwl=lwl,
                                          bwl=bwl, Tc=Tc,
                                          # hull_coe=None,
                                          keel_coe=(0.460, 0., - 0.2167),
                                          rudder_coe=(0.0536, 0., -0.11495),
                                          keel_sweep_back_angle=20.)

spkv09_swept_minus20 = \
    sideforce_production_keuning_verwerft(boatspeed=4.,
                                          heel_angle=0.,
                                          leeway_angle=4.,
                                          rudder_angle=0.,
                                          keel_projected_area=0.025240,
                                          keel_span=0.42 - Tc - bulb_diameter,
                                          rudder_projected_area=0.012402,
                                          rudder_span=0.248,
                                          lwl=lwl, bwl=bwl,
                                          Tc=Tc,
                                          # hull_coe=None,
                                          keel_coe=(0.460, 0., -0.2167),
                                          rudder_coe=(0.0536, 0., -0.11495),
                                          keel_sweep_back_angle=-20.)


def test_known_value_sideforce_production_keuning_verwerft():
    r"""No know data."""
    assert True


def test_negative_bwl_sideforce_production_keuning_verwerft():
    r"""Negative lwl."""
    with pytest.raises(ValueError):
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=0.,
                                              leeway_angle=4.,
                                              rudder_angle=0.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=-lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., - 0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))


def test_sign_and_symmetry_sideforce_production_keuning_verwerft():
    r"""Test signs and symmetry."""
    force = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=0.,
                                              leeway_angle=4.,
                                              rudder_angle=0.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))
    force_backwards = \
        sideforce_production_keuning_verwerft(boatspeed=-4.,
                                              heel_angle=0.,
                                              leeway_angle=4.,
                                              rudder_angle=0.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., - 0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))
    force_nose_down = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=0.,
                                              leeway_angle=-4.,
                                              rudder_angle=0.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., - 0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))

    assert force[0].fx + force[1].fx + force[2].fx == -(force_backwards[0].fx + force_backwards[1].fx + force_backwards[2].fx)
    assert force[0].fy + force[1].fy + force[2].fy == -(force_nose_down[0].fy + force_nose_down[1].fy + force_nose_down[2].fy)
    assert force[0].fx < 0
    assert force[1].fx < 0
    # print force[2].force[0]

    # hull drag is 0 in the model -> all the induced drag is on keel + rudder
    assert force[2].fx == 0


def test_heel_angle_resistance_sideforce_production_keuning_verwerft():
    r"""Test x force component with opposite heel angles."""
    force = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=20.,
                                              leeway_angle=4.,
                                              rudder_angle=3.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))
    force_heeled_to_windward_negative_leeway = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=-20.,
                                              leeway_angle=-4.,
                                              rudder_angle=-3.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248, lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., - 0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))
    assert force[0].fx + force[1].fx + force[2].fx == (force_heeled_to_windward_negative_leeway[0].fx +
                                                    force_heeled_to_windward_negative_leeway[1].fx +
                                                    force_heeled_to_windward_negative_leeway[2].fx)


def test_heel_angle_sideforce_sideforce_production_keuning_verwerft():
    r"""Test y force component with opposite heel angles."""
    force = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=20.,
                                              leeway_angle=4.,
                                              rudder_angle=3.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))
    force_heeled_to_windward_negative_leeway = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=-20.,
                                              leeway_angle=-4.,
                                              rudder_angle=-3.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., - 0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))
    assert force[0].fy + force[1].fy + force[2].fy == -(force_heeled_to_windward_negative_leeway[0].fy +
                                                     force_heeled_to_windward_negative_leeway[1].fy +
                                                     force_heeled_to_windward_negative_leeway[2].fy)


def test_rudder_angle_sideforce_production_keuning_verwerft():
    r"""Opposite rudder angles."""
    force = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=0.,
                                              leeway_angle=0.,
                                              rudder_angle=5.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))
    force_negative_rudder = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=0.,
                                              leeway_angle=0.,
                                              rudder_angle=-5.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495))
    # Reminder 0 -> keel, 1 -> rudder, 2 -> hull
    assert force[1].fx == force_negative_rudder[1].fx
    assert force[1].fy == -force_negative_rudder[1].fy


def test_sweep_symmetry_sideforce_production_keuning_verwerft():
    r"""Opposite keel sweep angles."""
    force = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=20.,
                                              leeway_angle=5.,
                                              rudder_angle=5.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495),
                                              keel_sweep_back_angle=20.)
    force_swept_minus = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=20.,
                                              leeway_angle=5.,
                                              rudder_angle=5.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., - 0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495),
                                              keel_sweep_back_angle=-20.)
    assert force[0].fx + force[1].fx + force[2].fx == (force_swept_minus[0].fx +
                                                    force_swept_minus[1].fx +
                                                    force_swept_minus[2].fx)

    assert force[0].fy + force[1].fy + force[2].fy == (force_swept_minus[0].fy +
                                                    force_swept_minus[1].fy +
                                                    force_swept_minus[2].fy)


def test_coe_position_sideforce_production_keuning_verwerft():
    r"""Centre of effort position tests."""
    sof_1 = SystemOfForces()
    sof_2 = SystemOfForces()

    force = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=20.,
                                              leeway_angle=5.,
                                              rudder_angle=5.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495),
                                              keel_sweep_back_angle=20.)
    for f in force:
        sof_1.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))

    assert sof_1.moment[2] < 0

    # assert force.position[0] > 0.0536
    # assert force.position[0] < 0.460
    # assert force.position[1] < 0.
    # assert force.position[2] > - 0.2167

    force = \
        sideforce_production_keuning_verwerft(boatspeed=4.,
                                              heel_angle=-20.,
                                              leeway_angle=5.,
                                              rudder_angle=5.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495),
                                              keel_sweep_back_angle=20.)
    for f in force:
        sof_2.add_force(Force((f.fx, f.fy, f.fz), (f.px, f.py, f.pz)))

    # assert force.position[0] > 0.0536
    # assert force.position[0] < 0.460
    # assert force.position[1] > 0.
    # assert force.position[2] > - 0.2167

    assert sof_2.moment[2] < 0


def test_zero_boatspeed_sideforce_production_keuning_verwerft():
    r"""The boatspeed is zero."""
    forces = \
        sideforce_production_keuning_verwerft(boatspeed=0.,
                                              heel_angle=-20.,
                                              leeway_angle=5.,
                                              rudder_angle=5.,
                                              keel_projected_area=0.025240,
                                              keel_span=0.42 - Tc - bulb_diameter,
                                              rudder_projected_area=0.012402,
                                              rudder_span=0.248,
                                              lwl=lwl,
                                              bwl=bwl,
                                              Tc=Tc,
                                              # hull_coe=None,
                                              keel_coe=(0.460, 0., -0.2167),
                                              rudder_coe=(0.0536, 0., -0.11495),
                                              keel_sweep_back_angle=20.)
    assert len(forces) == 3
    for f in forces:
        assert f.fx == 0.
        assert f.fy == 0.
        assert f.fz == 0.
