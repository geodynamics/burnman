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

    :param unit_cell_v: Unit cell volumes [A^3/unit cell].
    :type unit_cell_v: float

    :param z: Number of formula units per unit cell.
    :type z: float

    :returns: Volume [m^3/mol]
    :rtype: float
    """
    V = unit_cell_v * constants.Avogadro / 1.0e30 / z
    return V


def cell_parameters_to_vectors(cell_parameters, frame_convention):
    """
    Converts cell parameters to unit cell vectors.

    :param cell_parameters: An array containing the three lengths of the
        unit cell vectors [m], and the three angles [degrees].
        The first angle (:math:`\\alpha`) corresponds to the angle between the
        second and the third cell vectors, the second (:math:`\\beta`) to the
        angle between the first and third cell vectors, and the third
        (:math:`\\gamma`) to the angle between the first and second vectors.
    :type cell_parameters: numpy.array (1D)

    :param frame_convention: A list (c) defining the reference frame
        convention.  This function dictates that the c[0]th cell vector
        is colinear with the c[0]th axis, the c[1]th cell vector is
        perpendicular to the c[2]th axis,
        and the c[2]th cell vector is defined to give a right-handed
        coordinate system.
        In common crystallographic shorthand, x[c[0]] // a[c[0]],
        x[c[2]] // a[c[2]]^* (i.e. perpendicular to a[c[0]] and a[c[1]]).

    :type frame_convention: list of three integers

    :returns: The three vectors defining the parallelopiped cell [m].

    :rtype: numpy.array (2D)
    """
    lengths = cell_parameters[:3]
    angles_deg = cell_parameters[3:]
    angles = np.radians(angles_deg)

    c = frame_convention

    n2 = (np.cos(angles[c[0]]) - np.cos(angles[c[2]]) * np.cos(angles[c[1]])) / np.sin(
        angles[c[2]]
    )
    MT = np.zeros((3, 3))
    MT[c[0], c[0]] = lengths[c[0]]
    MT[c[1], c[0]] = lengths[c[1]] * np.cos(angles[c[2]])
    MT[c[1], c[1]] = lengths[c[1]] * np.sin(angles[c[2]])
    MT[c[2], c[0]] = lengths[c[2]] * np.cos(angles[c[1]])
    MT[c[2], c[1]] = lengths[c[2]] * n2
    MT[c[2], c[2]] = lengths[c[2]] * np.sqrt(np.sin(angles[c[1]]) ** 2 - n2**2)

    return MT


def cell_vectors_to_parameters(vectors, frame_convention):
    """
    Converts unit cell vectors to cell parameters.

    :param vectors: The three vectors defining the parallelopiped cell [m].
    :type vectors: numpy.array (2D)

    :param frame_convention: A list (c) defining the reference frame
        convention. This function dictates that the c[0]th cell vector
        is colinear with the c[0]th axis, the c[1]th cell vector is
        perpendicular to the c[2]th axis,
        and the c[2]th cell vector is defined to give a right-handed
        coordinate system.
        In common crystallographic shorthand, x[c[0]] // a[c[0]],
        x[c[2]] // a[c[2]]^* (i.e. perpendicular to a[c[0]] and a[c[1]]).

    :type frame_convention: list of three integers

    :returns: An array containing the three lengths of the unit cell vectors [m],
        and the three angles [degrees].
        The first angle (:math:`\\alpha`) corresponds to the angle between the
        second and the third cell vectors, the second (:math:`\\beta`) to the
        angle between the first and third cell vectors, and the third
        (:math:`\\gamma`) to the angle between the first and second vectors.
    :rtype: numpy.array (1D)
    """

    c = frame_convention
    assert np.abs(vectors[c[0], c[1]]) < np.finfo(float).eps
    assert np.abs(vectors[c[0], c[2]]) < np.finfo(float).eps
    assert np.abs(vectors[c[1], c[2]]) < np.finfo(float).eps

    lengths = np.empty(3)
    angles = np.empty(3)

    lengths[c[0]] = vectors[c[0], c[0]]
    lengths[c[1]] = np.sqrt(
        np.power(vectors[c[1], c[0]], 2.0) + np.power(vectors[c[1], c[1]], 2.0)
    )
    lengths[c[2]] = np.sqrt(
        np.power(vectors[c[2], c[0]], 2.0)
        + np.power(vectors[c[2], c[1]], 2.0)
        + np.power(vectors[c[2], c[2]], 2.0)
    )

    angles[c[2]] = np.arccos(vectors[c[1], c[0]] / lengths[c[1]])
    angles[c[1]] = np.arccos(vectors[c[2], c[0]] / lengths[c[2]])
    angles[c[0]] = np.arccos(
        vectors[c[2], c[1]] / lengths[c[2]] * np.sin(angles[c[2]])
        + np.cos(angles[c[2]]) * np.cos(angles[c[1]])
    )

    angles_deg = np.degrees(angles)
    return np.array([*lengths, *angles_deg])
