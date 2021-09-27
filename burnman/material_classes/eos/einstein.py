# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from .. import constants

"""
Functions for the Einstein model of a solid.
"""

eps = np.finfo(float).eps


def thermal_energy(T, einstein_T, n):
    """
    calculate the thermal energy of a substance.  Takes the temperature,
    the Einstein temperature, and n, the number of atoms per molecule.
    Returns thermal energy in J/mol
    """
    if T <= eps:
        # zero point energy
        return 3. * n * constants.gas_constant * einstein_T * 0.5
    x = einstein_T / T
    E_th = 3. * n * constants.gas_constant * einstein_T * \
        (0.5 + 1. / (np.exp(x) - 1.0))  # include the zero point energy
    return E_th


def molar_heat_capacity_v(T, einstein_T, n):
    """
    Heat capacity at constant volume.  In J/K/mol
    """
    if T <= eps:
        return 0.
    x = einstein_T / T
    C_v = 3.0 * n * constants.gas_constant * \
        (x * x * np.exp(x) / np.power(np.exp(x) - 1.0, 2.0))
    return C_v
