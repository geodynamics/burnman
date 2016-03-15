# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import numpy as np
from scipy.linalg import lu

from .processchemistry import *
from . import constants
from . import solidsolution

# This module computes chemical potentials (partial molar gibbs free
# energies) for an assemblage based on the Gibbs free energies and
# compositions of the individual phases.

# It can also calculate fugacities based on the gibbs free energies of the
# endmembers corresponding to chemical components.

# It can also calculate fugacities relative to other bulk compositions


def chemical_potentials(assemblage, component_formulae):
    """
    The compositional space of the components does not have to be a
    superset of the compositional space of the assemblage. Nor do they have to
    compose an orthogonal basis.

    The components must each be described by a linear mineral combination

    The mineral compositions must be linearly independent

        Parameters
        ----------
        assemblage : list of classes
            List of material classes
            set_method and set_state should already have been used
            the composition of the solid solutions should also have been set

        component_formulae : list of dictionaries
            List of chemical component formula dictionaries
            No restriction on length

        Returns
        -------
        component_potentials : array of floats
            Array of chemical potentials of components

    """
    # Split solid solutions into their respective endmembers
    # Find the chemical potentials of all the endmembers
    endmember_list = []
    endmember_potentials = []
    for mineral in assemblage:
        if isinstance(mineral, solidsolution.SolidSolution):
            for member in mineral.endmembers:
                endmember_list.append(member[0])
            for potential in mineral.partial_gibbs:
                endmember_potentials.append(potential)
        else:
            endmember_list.append(mineral)
            endmember_potentials.append(mineral.gibbs)

    # Make an array of all the endmember formulae
    endmember_formulae = [endmember.params['formula']
                          for endmember in endmember_list]
    endmember_compositions, elements = compositional_array(endmember_formulae)

    pl, u = lu(endmember_compositions, permute_l=True)
    assert(min(np.dot(np.square(u), np.ones(shape=len(elements)))) > 1e-6), \
        'Endmember compositions do not form an independent set of basis vectors'

    # Make an array of component formulae with elements in the same order as
    # the endmember array
    component_compositions = ordered_compositional_array(
        component_formulae, elements)

    p = np.linalg.lstsq(endmember_compositions.T, component_compositions.T)
    for idx, error in enumerate(p[1]):
        assert (error < 1e-10), \
            'Component %d not defined by prescribed assemblage' % (idx + 1)

    # Create an array of endmember proportions which sum to each component
    # composition
    endmember_proportions = np.around(p[0], 10).T

    # Calculate the chemical potential of each component
    component_potentials = np.dot(endmember_proportions, endmember_potentials)
    return component_potentials


def fugacity(standard_material, assemblage):
    """
        Parameters
        ----------
        standard_material: class
            Material class
            set_method and set_state should already have been used
            material must have a formula as a dictionary parameter

        assemblage: list of classes
            List of material classes
            set_method and set_state should already have been used

        Returns
        -------
        fugacity : float
            Value of the fugacity of the component with respect to
            the standard material

    """
    component_formula = standard_material.params['formula']
    chemical_potential = chemical_potentials(
        assemblage, [component_formula])[0]

    fugacity = np.exp((chemical_potential - standard_material.gibbs) / (
        constants.gas_constant * assemblage[0].temperature))
    return fugacity


def relative_fugacity(standard_material, assemblage, reference_assemblage):
    """
        Parameters
        ----------
        standard_material: class
            Material class
            set_method and set_state should already have been used
            material must have a formula as a dictionary parameter

        assemblage: list of classes
            List of material classes
            set_method and set_state should already have been used

        reference_assemblage: list of classes
            List of material classes
            set_method and set_state should already have been used

        Returns
        -------
        relative_fugacity : float
            Value of the fugacity of the component in the assemblage
            with respect to the reference_assemblage

    """
    component_formula = standard_material.params['formula']
    chemical_potential = chemical_potentials(
        assemblage, [component_formula])[0]
    reference_chemical_potential = chemical_potentials(
        reference_assemblage, [component_formula])[0]

    relative_fugacity = np.exp((chemical_potential - reference_chemical_potential) / (
        constants.gas_constant * assemblage[0].temperature))
    return relative_fugacity
