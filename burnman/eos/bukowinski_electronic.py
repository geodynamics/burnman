# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import numpy as np

"""
Functions for the Bukowinski (1977) model for the electronic component
of the Helmholtz energy as used by Stixrude an Lithgow-Bertelloni (2024).
"""


def helmholtz(temperature, volume, T_0, V_0, bel_0, gel):
    return (
        -0.5
        * bel_0
        * np.power(volume / V_0, gel)
        * (temperature * temperature - T_0 * T_0)
    )


def pressure(temperature, volume, T_0, V_0, bel_0, gel):
    """
    P = -dF/dV
    """
    return (
        0.5
        * gel
        * bel_0
        * np.power(volume / V_0, gel)
        * (temperature * temperature - T_0 * T_0)
        / volume
    )


def KToverV(temperature, volume, T_0, V_0, bel_0, gel):
    """
    KT = -V dP/dV
    """
    return -(gel - 1.0) * pressure(temperature, volume, T_0, V_0, bel_0, gel) / volume


def entropy(temperature, volume, V_0, bel_0, gel):
    """
    S = -dF/dT
    """
    return bel_0 * np.power(volume / V_0, gel) * temperature


def CVoverT(volume, V_0, bel_0, gel):
    """
    CV = T dS/dT
    """
    return bel_0 * np.power(volume / V_0, gel)


def aKT(temperature, volume, V_0, bel_0, gel):
    """
    aKT = dP/dT
    """
    return gel * bel_0 * np.power(volume / V_0, gel) * temperature / volume
