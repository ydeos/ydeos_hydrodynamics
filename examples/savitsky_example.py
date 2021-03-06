#!/usr/bin/env python
# coding: utf-8

r"""Savitsky example.

Values from p204 example in Principles of Yacht Design

"""

from ydeos_hydrodynamics.savitsky import savitsky_report

savitsky_report(boatspeed=20.6,
                m=9750.,
                LCG=3.55,
                VCG=1.145,
                b=2.74,
                epsilon=3.,
                beta=17.9,
                f=0.49)
