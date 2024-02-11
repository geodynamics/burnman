# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import numpy as np
from scipy.optimize import brentq
from ..utils.misc import read_table, lookup_and_interpolate
from ..utils.math import bracket


def brown_shankland(depths):
    """
    Geotherm from :cite:`Brown1981`. NOTE: Valid only above 270 km

    :param depths: The list of depths at which to evaluate the geotherm.
        :math:`[m]`
    :type depths: list of floats

    :returns: The list of temperatures for each of the pressures. :math:`[K]`
    :rtype: list of floats
    """

    assert min(depths) >= min(table_brown_depth)
    assert max(depths) <= max(table_brown_depth)
    temperature = np.empty_like(depths)
    for i, depth in enumerate(depths):
        temperature[i] = lookup_and_interpolate(
            table_brown_depth, table_brown_temperature, depth
        )
    return temperature


def anderson(depths):
    """
    Geotherm from :cite:`anderson1982earth`.

    :param depths: The list of depths at which to evaluate the geotherm.
        :math:`[m]`
    :type depths: list of floats

    :returns: The list of temperatures for each of the pressures. :math:`[K]`
    :rtype: list of floats
    """
    assert min(depths) >= min(table_anderson_depth)
    assert max(depths) <= max(table_anderson_depth)
    temperature = np.empty_like(depths)
    for i, depth in enumerate(depths):
        temperature[i] = lookup_and_interpolate(
            table_anderson_depth, table_anderson_temperature, depth
        )
    return temperature


def stacey_continental(depths):
    """
    Continental geotherm from :cite:`stacey1977`.

    :param depths: The list of depths at which to evaluate the geotherm.
        :math:`[m]`
    :type depths: list of floats

    :returns: The list of temperatures for each of the pressures. :math:`[K]`
    :rtype: list of floats
    """
    assert min(depths) >= min(table_stacey_c_depth)
    assert max(depths) <= max(table_stacey_c_depth)
    temperature = np.empty_like(depths)
    for i, depth in enumerate(depths):
        temperature[i] = lookup_and_interpolate(
            table_stacey_c_depth, table_stacey_c_temperature, depth
        )
    return temperature


def stacey_oceanic(depths):
    """
    Oceanic geotherm from :cite:`stacey1977`.

    :param depths: The list of depths at which to evaluate the geotherm.
        :math:`[m]`
    :type depths: list of floats

    :returns: The list of temperatures for each of the pressures. :math:`[K]`
    :rtype: list of floats
    """
    assert min(depths) >= min(table_stacey_o_depth)
    assert max(depths) <= max(table_stacey_o_depth)
    temperature = np.empty_like(depths)
    for i, depth in enumerate(depths):
        temperature[i] = lookup_and_interpolate(
            table_stacey_o_depth, table_stacey_o_temperature, depth
        )
    return temperature


def adiabatic(pressures, T0, rock):
    """
    This calculates a geotherm based on an anchor temperature and a rock,
    assuming that the rock's temperature follows an adiabatic gradient with
    pressure. A good first guess is provided by integrating:

    .. math::
        \\frac{\\partial T}{\\partial P} = \\frac{ \\gamma  T}{ K_s }

    where :math:`\\gamma` is the Grueneisen parameter and :math:`K_s` is
    the adiabatic bulk modulus.

    :param pressures: The list of pressures in :math:`[Pa]` at which
        to evaluate the geotherm.
    :type pressures: list of floats

    :param T0: An anchor temperature, corresponding to the temperature of the first
        pressure in the list. :math:`[K]`
    :type T0: float

    :param rock: Composite for which we compute the adiabat.  From this material we
        must compute average Grueneisen parameters and adiabatic bulk moduli
        for each pressure/temperature.
    :type rock: :class:`burnman.composite`

    :returns: The list of temperatures for each pressure. :math:`[K]`
    :rtype: numpy.array of floats
    """

    rock.set_state(pressures[0], T0)
    S0 = rock.S

    delta_S = lambda T, P, rock, S0: S0 - rock.evaluate(["S"], [P], [T])[0][0]

    temperatures = np.empty_like(pressures)
    temperatures[0] = T0
    for i in range(1, len(pressures)):
        args = (pressures[i], rock, S0)
        sol = bracket(
            fn=delta_S,
            x0=(
                temperatures[i - 1]
                + (rock.gr * temperatures[i - 1] / rock.K_S)
                * (pressures[i] - pressures[i - 1])
            ),
            dx=1.0,
            args=args,
        )
        temperatures[i] = brentq(delta_S, sol[0], sol[1], args=args)

    return temperatures


table_brown = read_table("input_geotherm/brown_81.txt")
table_brown_depth = np.array(table_brown)[:, 0]
table_brown_temperature = np.array(table_brown)[:, 1]

table_anderson = read_table("input_geotherm/anderson_82.txt")
table_anderson_depth = np.array(table_anderson)[:, 0]
table_anderson_temperature = np.array(table_anderson)[:, 1]

table_stacey_c = read_table("input_geotherm/Stacey_1977_continents.txt")
table_stacey_c_depth = np.array(table_stacey_c)[:, 0]
table_stacey_c_temperature = np.array(table_stacey_c)[:, 1]

table_stacey_o = read_table("input_geotherm/Stacey_1977_oceans.txt")
table_stacey_o_depth = np.array(table_stacey_o)[:, 0]
table_stacey_o_temperature = np.array(table_stacey_o)[:, 1]
