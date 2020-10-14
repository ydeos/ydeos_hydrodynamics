# coding: utf-8

r"""Form factors."""

from math import sqrt


def k_watanabe(lwl: float, bwl: float, Tc: float, block_coefficient: float) -> float:
    r"""Watanabe estimate of the k factor.

    As defined in Resistance and Wake Prediction
    for Early Stage Ship Design by Brian Johnson
    lwl, bwl, Tc : [m]
    Reminder : Cv = (1 + k) * Cf

    """
    return -0.095 + 25.6 * (block_coefficient / ((lwl/bwl)**2 * sqrt(bwl/Tc)))


def hoerner(thickness_to_chord: float) -> float:
    r"""Hoerner form factor for a foil.

    From thickness to chord ratio of 2D foil section

    """
    return 1 + 2 * thickness_to_chord + 60 * thickness_to_chord ** 4


def hoerner_streamlined_body(thickness_to_chord: float) -> float:
    r"""Hoerner form factor for streamlined bodies at Rn > 1e5."""
    return 1 + 1.5 * thickness_to_chord ** 1.5 + 7 * thickness_to_chord ** 3
