# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import numpy as np

# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    import os
    if 'NUMBA_DISABLE_JIT' in os.environ and int(os.environ['NUMBA_DISABLE_JIT']) == 1:
        raise ImportError("NOOOO!")
    from numba import jit
except ImportError:
    def jit(fn):
        return fn

import scipy.integrate as integrate

from .. import constants

"""
Functions for the Debye model.  Note that this is not Mie-Grueneisen-Debye,
just Debye, so is pretty limited.  Combine this with Mie-Grueneisen and
Birch-Murnaghan to get a full EOS
"""


chebyshev_representation = np.array(
    [2.707737068327440945 / 2.0, 0.340068135211091751, -0.12945150184440869e-01,
     0.7963755380173816e-03, -
     0.546360009590824e-04, 0.39243019598805e-05,
     -0.2894032823539e-06, 0.217317613962e-07, -
     0.16542099950e-08,
     0.1272796189e-09, -
     0.987963460e-11, 0.7725074e-12, -
     0.607797e-13,
     0.48076e-14, -0.3820e-15, 0.305e-16, -0.24e-17])


@jit
def _chebval(x, c):
    """
    Evaluate a Chebyshev series at points x.
    This is just a lightly modified copy/paste job from the numpy
    implementation of the same function, copied over here to put a
    jit wrapper around it.
    """
    if len(c) == 1:
        c0 = c[0]
        c1 = 0
    elif len(c) == 2:
        c0 = c[0]
        c1 = c[1]
    else:
        x2 = 2 * x
        c0 = c[-2]
        c1 = c[-1]
        for i in range(3, len(c) + 1):
            tmp = c0
            c0 = c[-i] - c1
            c1 = tmp + c1 * x2
    return c0 + c1 * x


def debye_fn(x):
    """
    Evaluate the Debye function.  Takes the parameter
    xi = Debye_T/T
    """
    sol = integrate.quad(
        lambda xi: (xi * xi * xi) / (np.exp(xi) - 1.), 0.0, x)  # EQ B3
    return 3. * sol[0] / pow(x, 3.)


eps = np.finfo(np.float).eps
sqrt_eps = np.sqrt(np.finfo(np.float).eps)
log_eps = np.log(np.finfo(np.float).eps)


@jit
def debye_fn_cheb(x):
    """
    Evaluate the Debye function using a Chebyshev series expansion coupled with
    asymptotic solutions of the function.  Shamelessly adapted from the GSL implementation
    of the same function (Itself adapted from Collected Algorithms from ACM).
    Should give the same result as debye_fn(x) to near machine-precision.
    """
    val_infinity = 19.4818182068004875
    xcut = -log_eps

    assert(x > 0.0)  # check for invalid x

    if x < 2.0 * np.sqrt(2.0) * sqrt_eps:
        return 1.0 - 3.0 * x / 8.0 + x * x / 20.0
    elif x <= 4.0:
        t = x * x / 8.0 - 1.0
        c = _chebval(t, chebyshev_representation)
        return c - 0.375 * x
    elif x < -(np.log(2.0) + log_eps):
        nexp = int(np.floor(xcut / x))
        ex = np.exp(-x)
        xk = nexp * x
        rk = nexp
        sum = 0.0
        for i in range(nexp, 0, -1):
            xk_inv = 1.0 / xk
            sum *= ex
            sum += (((6.0 * xk_inv + 6.0) * xk_inv + 3.0) * xk_inv + 1.0) / rk
            rk -= 1.0
            xk -= x
        return val_infinity / (x * x * x) - 3.0 * sum * ex
    elif x < xcut:
        x3 = x * x * x
        sum = 6.0 + 6.0 * x + 3.0 * x * x + x3
        return (val_infinity - 3.0 * sum * np.exp(-x)) / x3
    else:
        return ((val_infinity / x) / x) / x


@jit
def thermal_energy(T, debye_T, n):
    """
    calculate the thermal energy of a substance.  Takes the temperature,
    the Debye temperature, and n, the number of atoms per molecule.
    Returns thermal energy in J/mol
    """
    if T <= eps:
        return 0.
    E_th = 3. * n * constants.gas_constant * T * debye_fn_cheb(debye_T / T)
    return E_th


@jit
def heat_capacity_v(T, debye_T, n):
    """
    Heat capacity at constant volume.  In J/K/mol
    """
    if T <= eps:
        return 0.
    x = debye_T / T
    C_v = 3.0 * n * constants.gas_constant * \
        (4.0 * debye_fn_cheb(x) - 3.0 * x / (np.exp(x) - 1.0))
    return C_v


@jit
def helmholtz_free_energy(T, debye_T, n):
    """
    Helmholtz free energy of lattice vibrations in the Debye model.
    It is important to note that this does NOT include the zero
    point energy of vibration for the lattice.  As long as you are
    calculating relative differences in F, this should cancel anyways.
    In Joules.
    """
    if T <= eps:
        return 0.
    x = debye_T / T
    F = n * constants.gas_constant * T * \
        (3.0 * np.log(1.0 - np.exp(-x)) - debye_fn_cheb(x))
    return F


def entropy(T, debye_T, n):
    """
    Entropy due to lattice vibrations in the Debye model [J/K]
    """
    if T <= eps:
        return 0.
    x = debye_T / T
    S = n * constants.gas_constant * \
        (4. * debye_fn_cheb(x) - 3. * np.log(1.0 - np.exp(-x)))
    return S
