# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import numpy as np
import scipy.integrate as integrate
from . import tools
from . import seismic


def brown_shankland(pressure):
    """
    Geotherm from :cite:`Brown1981`. NOTE: Valid only above 270 km

    Parameters
    ----------
    pressure : list of floats
        The list of pressures at which to evaluate the geotherm. :math:`[Pa]`

    Returns
    -------
    temperature : list of floats
        The list of temperatures for each of the pressures. :math:`[K]`
    """
    temperature = np.empty_like(pressure)
    for i in range(len(pressure)):
        depth = seismic.prem_model.depth(pressure[i])
        if depth < min(table_brown_depth):
            raise ValueError(
                "depth smaller than range Brown & Shankland, 1981")
        temperature[i] = tools.lookup_and_interpolate(
            table_brown_depth, table_brown_temperature, depth)
    return temperature


def anderson(pressure):
    """
    Geotherm from :cite:`anderson1982earth`.

    Parameters
    ----------
    pressure : list of floats
        The list of pressures at which to evaluate the geotherm. :math:`[Pa]`

    Returns
    -------
    temperature : list of floats
        The list of temperatures for each of the pressures. :math:`[K]`
    """
    temperature = np.empty_like(pressure)
    for i in range(len(pressure)):
        depth = seismic.prem_model.depth(pressure[i])
        temperature[i] = tools.lookup_and_interpolate(
            table_anderson_depth, table_anderson_temperature, depth)
    return temperature


def adiabatic(pressures, T0, rock):
    """
    This calculates a geotherm based on an anchor temperature and a rock,
    assuming that the rock's temperature follows an adiabatic gradient with
    pressure. This amounts to integrating:

    .. math::
        \\frac{\partial T}{\partial P} = \\frac{ \\gamma  T}{ K_s }

    where :math:`\\gamma` is the Grueneisen parameter and :math:`K_s` is
    the adiabatic bulk modulus.

    Parameters
    ----------

    pressures : list of floats
        The list of pressures in :math:`[Pa]` at which to evaluate the geotherm.

    T0 : float
        An anchor temperature, corresponding to the temperature of the first
        pressure in the list. :math:`[K]`

    rock : :class:`burnman.composite`
        Material for which we compute the adiabat.  From this material we
        must compute average Grueneisen parameters and adiabatic bulk moduli
        for each pressure/temperature.

    Returns
    -------

    temperature: list of floats
        The list of temperatures for each pressure. :math:`[K]`
    """
    temperatures = integrate.odeint(
        lambda t, p: dTdP(t, p, rock), T0, pressures)
    return temperatures.ravel()


def dTdP(temperature, pressure, rock):
    """
    ODE to integrate temperature with depth for a composite material
    Assumes that the minerals exist at a common pressure (Reuss bound, should be good for
    slow deformations at high temperature), as well as an adiabatic process.  This
    corresponds to conservation of enthalpy.
    First consider compression of the composite to a new pressure P+dP.  They all heat up
    different amounts dT[i], according to their thermoelastic parameters.  Then allow them
    to equilibrate to a constant temperature dT, conserving heat within the composite.
    This works out to the formula:

    .. math::
        dT/dP = T*\\frac{\Sigma_i(X[i]*C_{p}[i]*\gamma[i]/K[i])}{\Sigma(X[i]*C_{p}[i])}

    Where :math:`X[i]` is the molar fraction of phase :math:`i`, :math:`C_p` is the specific heat at constant pressure,
    :math:`\gamma` is the Gruneisen parameter and :math:`K` is the bulk modulus.
    This function is called by :func:`burnman.geotherm.adiabatic`, and in general
    it will not be too useful in other contexts.

    Parameters
    ----------

    pressure : float
        The pressure at which to evaluate dT/dP. :math:`[Pa]`

    temperature : float
        The temperature at which to evaluate dT/dP. :math:`[K]`

    rock : :class:`burnman.composite`
        Material for which we compute dT/dP.

    Returns
    -------
        dT/dP : float
          Adiabatic temperature gradient for the composite at a given temperature and pressure. :math:`[K/Pa]`
    """
    top = 0
    bottom = 0
    rock.set_state(pressure, temperature)
    (minerals, fractions) = rock.unroll()
    for (mineral, fraction) in zip(minerals, fractions):
        gr = mineral.grueneisen_parameter
        K_s = mineral.adiabatic_bulk_modulus
        C_p = mineral.heat_capacity_p

        top += fraction * gr * C_p / K_s
        bottom += fraction * C_p

    return temperature * top / bottom


table_brown = tools.read_table("input_geotherm/brown_81.txt")
table_brown_depth = np.array(table_brown)[:, 0]
table_brown_temperature = np.array(table_brown)[:, 1]

table_anderson = tools.read_table("input_geotherm/anderson_82.txt")
table_anderson_depth = np.array(table_anderson)[:, 0]
table_anderson_temperature = np.array(table_anderson)[:, 1]
