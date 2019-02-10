# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
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
