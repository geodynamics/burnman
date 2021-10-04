# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from .. import constants


def molar_volume_from_unit_cell_volume(unit_cell_v, z):
    """
    Converts a unit cell volume from Angstroms^3 per unitcell,
    to m^3/mol.

    Parameters
    ----------
    unit_cell_v : float
        Unit cell volumes [A^3/unit cell]

    z : float
        Number of formula units per unit cell


    Returns
    -------
    V : float
        Volume [m^3/mol]
    """
    V = unit_cell_v * constants.Avogadro / 1.e30 / z
    return V


def cell_parameters_to_vectors(cell_parameters):
    """
    Converts cell parameters to unit cell vectors.

    Parameters
    ----------
    cell_parameters : 1D numpy array
        An array containing the three lengths of the unit cell vectors [m],
        and the three angles [degrees].
        The first angle (:math:`\\alpha`) corresponds to the angle between the
        second and the third cell vectors, the second (:math:`\\beta`) to the
        angle between the first and third cell vectors, and the third
        (:math:`\\gamma`) to the angle between the first and second vectors.

    Returns
    -------
    M : 2D numpy array
        The three vectors defining the parallelopiped cell [m].
        This function assumes that the first cell vector is colinear with the
        x-axis, and the second is perpendicular to the z-axis, and the third is
        defined in a right-handed sense.
    """
    a, b, c, alpha_deg, beta_deg, gamma_deg = cell_parameters
    alpha = np.radians(alpha_deg)
    beta = np.radians(beta_deg)
    gamma = np.radians(gamma_deg)

    n2 = (np.cos(alpha)-np.cos(gamma)*np.cos(beta))/np.sin(gamma)
    M = np.array([[a, 0, 0],
                  [b*np.cos(gamma), b*np.sin(gamma), 0],
                  [c*np.cos(beta), c*n2, c*np.sqrt(np.sin(beta)**2-n2**2)]])
    return M


def cell_vectors_to_parameters(M):
    """
    Converts unit cell vectors to cell parameters.

    Parameters
    ----------
    M : 2D numpy array
        The three vectors defining the parallelopiped cell [m].
        This function assumes that the first cell vector is colinear with the
        x-axis, the second is perpendicular to the z-axis, and the third is
        defined in a right-handed sense.

    Returns
    -------
    cell_parameters : 1D numpy array
        An array containing the three lengths of the unit cell vectors [m],
        and the three angles [degrees].
        The first angle (:math:`\\alpha`) corresponds to the angle between the
        second and the third cell vectors, the second (:math:`\\beta`) to the
        angle between the first and third cell vectors, and the third
        (:math:`\\gamma`) to the angle between the first and second vectors.
    """

    assert M[0, 1] == 0
    assert M[0, 2] == 0
    assert M[1, 2] == 0

    a = M[0, 0]
    b = np.sqrt(np.power(M[1, 0], 2.) + np.power(M[1, 1], 2.))
    c = (np.sqrt(np.power(M[2, 0], 2.)
                 + np.power(M[2, 1], 2.)
                 + np.power(M[2, 2], 2.)))

    gamma = np.arccos(M[1, 0] / b)
    beta = np.arccos(M[2, 0] / c)
    alpha = np.arccos(M[2, 1]/c*np.sin(gamma) + np.cos(gamma)*np.cos(beta))

    gamma_deg = np.degrees(gamma)
    beta_deg = np.degrees(beta)
    alpha_deg = np.degrees(alpha)

    return np.array([a, b, c, alpha_deg, beta_deg, gamma_deg])
