# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import numpy as np
import scipy.integrate as integrate

from . import geotherm
from . import seismic
from . import averaging_schemes


def velocities_from_rock(rock, pressures, temperatures, averaging_scheme=averaging_schemes.VoigtReussHill()):
    """
    This function is deprecated. Use :func:`burnman.material.Material.evaluate` instead.

    A function that rolls several steps into one: given a rock and a list of
    pressures and temperatures, it calculates the elastic moduli of the
    individual phases using calculate_moduli(), averages them using
    average_moduli(), and calculates the seismic velocities using
    compute_velocities().

    Parameters
    ----------
    rock : :class:`burnman.Material`
         this is the rock for which you are calculating velocities

    pressures: list of float
        list of pressures you want to evaluate the rock at. :math:`[Pa]`

    temperatures: list of float
        list of temperatures you want to evaluate the rock at. :math:`[K]`


    averaging_scheme: :class:`burnman.averaging_schemes.AveragingScheme`
        Averaging scheme to use.

    Returns
    -------

    rho, V_p, V_s, V_phi, K, G : lists of floats
        Lists of density [kg/m^3], P-wave velocity [m/s], shear-wave velocity [m/s], bulk sound velocity [m/s], bulk modulus [Pa], and shear modulus [Pa] for each P,T point.

    """
    old_averaging_scheme = rock.averaging_scheme
    rock.set_averaging_scheme(averaging_scheme)
    rho, vp, vs, vphi, K, G = rock.evaluate(
        ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], pressures, temperatures)
    rock.set_averaging_scheme(old_averaging_scheme)
    return rho, vp, vs, vphi, K, G


def compare_l2(depth, calc, obs):
    """
    PUT IN TOOLS
    Computes the L2 norm for N profiles at a time (assumed to be linear between points).

    .. math:: math does not work yet...
       \sum_{i=1}^{\\infty} x_{i}

    :type depths: array of float
    :param depths: depths. :math:`[m]`
    :type calc: list of arrays of float
    :param calc: N arrays calculated values, e.g. [mat_vs,mat_vphi]
    :type obs: list of arrays of float
    :param obs: N arrays of values (observed or calculated) to compare to , e.g. [seis_vs, seis_vphi]

    :returns: array of L2 norms of length N
    :rtype: array of floats
    """
    err = []
    for l in range(len(calc)):
        err.append(l2(depth, calc[l], obs[l]))

    return err


def compare_chifactor(calc, obs):
    """
    PUT IN TOOLS
    Computes the chi factor for N profiles at a time. Assumes a 1% a priori uncertainty on the seismic model.


    :type calc: list of arrays of float
    :param calc: N arrays calculated values, e.g. [mat_vs,mat_vphi]
    :type obs: list of arrays of float
    :param obs: N arrays of values (observed or calculated) to compare to , e.g. [seis_vs, seis_vphi]

    :returns: error array of length N
    :rtype: array of floats
    """
    err = []
    for l in range(len(calc)):
        err.append(chi_factor(calc[l], obs[l]))

    return err


def l2(x, funca, funcb):
    """
    PUT IN TOOLS
    Computes the L2 norm for one profile(assumed to be linear between points).

    :type x: array of float
    :param x: depths :math:`[m]`.
    :type funca: list of arrays of float
    :param funca: array calculated values
    :type funcb: list of arrays of float
    :param funcb: array of values (observed or calculated) to compare to

    :returns: L2 norm
    :rtype: array of floats
    """
    diff = np.array(funca - funcb)
    diff = diff * diff
    return integrate.trapz(diff, x)


def nrmse(x, funca, funcb):
    """
    PUT IN TOOLS
    Normalized root mean square error for one profile
    :type x: array of float
    :param x: depths in m.
    :type funca: list of arrays of float
    :param funca: array calculated values
    :type funcb: list of arrays of float
    :param funcb: array of values (observed or calculated) to compare to

    :returns: RMS error
    :rtype: array of floats

    """
    diff = np.array(funca - funcb)
    diff = diff * diff
    rmse = np.sqrt(np.sum(diff) / x)
    nrmse = rmse / (np.max(funca) - np.min(funca))
    return nrmse


def chi_factor(calc, obs):
    """
    PUT IN TOOLS
    :math:`\\chi` factor for one profile assuming 1% uncertainty on the reference model (obs)
    :type calc: list of arrays of float
    :param calc: array calculated values
    :type obs: list of arrays of float
    :param obs: array of reference values to compare to

    :returns: :math:`\\chi` factor
    :rtype: array of floats

    """

    err = np.empty_like(calc)
    for i in range(len(calc)):
        err[i] = pow((calc[i] - obs[i]) / (0.01 * np.mean(obs)), 2.)

    err_tot = np.sum(err) / len(err)

    return err_tot
