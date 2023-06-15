# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from .. import constants

# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    import os

    if "NUMBA_DISABLE_JIT" in os.environ and int(os.environ["NUMBA_DISABLE_JIT"]) == 1:
        raise ImportError("NOOOO!")
    from numba import jit
except ImportError:

    def jit(fn):
        return fn


"""
Functions for the Einstein model of a solid.
"""

eps = np.finfo(float).eps


@jit(nopython=True)
def thermal_energy(T, einstein_T, n):
    """
    calculate the thermal energy of a substance.  Takes the temperature,
    the Einstein temperature, and n, the number of atoms per molecule.
    Returns thermal energy in J/mol
    """
    if T <= eps:
        return 0.0
    x = einstein_T / T
    E_th = 3.0 * n * constants.gas_constant * einstein_T * (1.0 / (np.exp(x) - 1.0))
    return E_th


@jit(nopython=True)
def molar_heat_capacity_v(T, einstein_T, n):
    """
    Heat capacity at constant volume.  In J/K/mol
    """
    if T <= eps:
        return 0.0
    x = einstein_T / T
    C_v = (
        3.0
        * n
        * constants.gas_constant
        * (x * x * np.exp(x) / np.power(np.exp(x) - 1.0, 2.0))
    )
    return C_v


@jit(nopython=True)
def helmholtz_free_energy(T, einstein_T, n):
    """
    Helmholtz free energy of lattice vibrations in the Einstein model [J].
    It is important to note that this does NOT include the zero
    point energy for the lattice.  As long as you are
    calculating relative differences in F, this should cancel anyway.
    """
    E = thermal_energy(T, einstein_T, n)
    S = entropy(T, einstein_T, n)
    return E - T * S


@jit(nopython=True)
def entropy(T, einstein_T, n):
    """
    Entropy due to lattice vibrations in the Einstein model [J/K]
    """
    if T <= eps:
        return 0.0
    x = einstein_T / T
    S = (
        3.0
        * n
        * constants.gas_constant
        * (-x * np.exp(-x) / (np.exp(-x) - 1.0) - np.log(1.0 - np.exp(-x)))
    )
    return S


@jit(nopython=True)
def dmolar_heat_capacity_v_dT(T, einstein_T, n):
    """
    First temperature derivative of the heat capacity at constant volume
    according to the Einstein model [J/K^2/mol].
    """
    if T <= eps:
        return 0.0

    x = einstein_T / T
    dCvdT = (
        3.0
        * n
        * constants.gas_constant
        * x
        * x
        * np.exp(x)
        * ((x - 2.0) * np.exp(x) + (x + 2.0))
        / (T * np.power(np.exp(x) - 1.0, 3.0))
    )

    return dCvdT
