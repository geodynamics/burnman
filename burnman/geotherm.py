# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import matplotlib.pyplot as pyplot
import scipy.integrate as integrate

import burnman.tools
import burnman.seismic


def brown_shankland(pressure):
    """
    Geotherm from Brown and Shankland (1981).

    Parameters
    ----------
    pressure : list of floats
        The list of pressures in [Pa] at which to evaluate the geotherm.

    Returns
    -------
    temperature : list of floats
        The list of temperatures in [K] for each of the pressures
    """
    temperature = np.empty_like(pressure)
    for i in range(len(pressure)):
      depth = burnman.seismic.prem_model.depth(pressure[i])
      temperature[i] = burnman.tools.lookup_and_interpolate(table_brown_depth, table_brown_temperature, depth)
    return temperature


def anderson(pressure):
    """
    Geotherm from Anderson (1982).

    Parameters
    ----------
    pressure : list of floats
        The list of pressures in [Pa] at which to evaluate the geotherm.

    Returns
    -------
    temperature : list of floats
        The list of temperatures in [K] for each of the pressures
    """
    temperature = np.empty_like(pressure)
    for i in range(len(pressure)):
      depth = burnman.seismic.prem_model.depth(pressure[i])
      temperature[i] = burnman.tools.lookup_and_interpolate(table_anderson_depth, table_anderson_temperature, depth)
    return temperature

def adiabatic(pressures, T0, rock):
    """
    This calculates a geotherm based on an anchor temperature and a rock,
    assuming that the rock's temperature follows an adiabatic gradient with
    pressure.  This amounts to integrating

    .. math::
        \\frac{\partial T}{\partial P} = \\frac{ \\gamma  T}{ K_s }

    where :math:`\\gamma` is the Grueneisen parameter and :math:`K_s` is
    the adiabatic bulk modulus.

    Parameters
    ----------

    pressures : list of floats
        The list of pressures in [Pa] at which to evaluate the geotherm.

    T0 : float
        An anchor temperature, corresponding to the temperature of the first
        pressure in the list.

    rock : :class:`burnman.composite`
        Material for which we compute the adiabat.  From this material we
        must compute average Grueneisen parameters and adiabatic bulk moduli
        for each pressure/temperature.

    Returns
    -------

    temperature: list of floats
        The list of temperatures in [K] for each of the pressures
    """
    temperatures = integrate.odeint(lambda t,p : dTdP(t,p,rock), T0, pressures)
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
    This works out to the formula: dT/dP = T*sum(frac[i]*Cp[i]*gr[i]/K[i])/sum(frac[i]*Cp[i])

    This function is called by :func:`burnman.geotherm.adiabatic`, and in general
    it will not be too useful in other contexts.

    Parameters
    ----------

    pressure : float
        The pressure at which to evaluate dT/dP.

    temperature : float
        The temperature at which to evaluate dT/dP.

    rock : :class:`burnman.composite`
        Material for which we compute dT/dP

    Returns
    -------
        dT/dP : float
          Adiabatic temperature gradient [K/Pa] for the composite at temperature, pressure.
    """
    top = 0
    bottom = 0
    rock.set_state(pressure, temperature)
    (fractions,minerals) = rock.unroll()
    for (fr,mineral) in zip(fractions,minerals):
        gr = mineral.grueneisen_parameter()
        K_s = mineral.adiabatic_bulk_modulus()
        C_p = mineral.heat_capacity_p()

        top += fr*gr*C_p/K_s
        bottom += fr*C_p

    return temperature*top/bottom


table_brown = burnman.tools.read_table("input_geotherm/brown_81.txt")
table_brown_depth = np.array(table_brown)[:,0]
table_brown_temperature = np.array(table_brown)[:,1]

table_anderson = burnman.tools.read_table("input_geotherm/anderson_82.txt")
table_anderson_depth = np.array(table_anderson)[:,0]
table_anderson_temperature = np.array(table_anderson)[:,1]

