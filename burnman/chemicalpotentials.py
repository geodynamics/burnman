# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import burnman
from burnman import minerals
from processchemistry import *
import numpy as np
from scipy.linalg import lu

R=8.3145
# This module computes chemical potentials (partial molar gibbs free energies) for an assemblage based on the Gibbs free energies and compositions of the individual phases.

# It can also calculate fugacities based on the gibbs free energies of the endmembers corresponding to chemical components.

# It can also calculate fugacities relative to other bulk compositions

def chemicalpotentials(assemblage, component_formulae):
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

        component_formulae : list of dictionaries
            List of chemical component formula dictionaries
            No restriction on length

        Returns
        -------
        component_potentials : array of floats
            Array of chemical potentials of components

    """
    mineral_gibbs=[mineral.gibbs for mineral in assemblage]
    mineral_formulae=[mineral.params['formula'] for mineral in assemblage]
    mineral_compositions, elements=compositional_array(mineral_formulae)

    pl, u = lu(mineral_compositions, permute_l=True)
    assert(min(np.dot(np.square(u), np.ones(shape=len(elements)))) > 1e-6), \
        'Mineral compositions do not form an independent set of basis vectors'

    component_compositions=ordered_compositional_array(component_formulae, elements)

    p=np.linalg.lstsq(mineral_compositions.T,component_compositions.T)
    for idx, error in enumerate(p[1]):
        assert (error < 1e-10), \
            'Component %d not defined by prescribed assemblage' % (idx+1)

    mineral_proportions=np.around(p[0],10).T

    component_potentials=np.dot(mineral_proportions, mineral_gibbs)

    return component_potentials

def fugacity(component_formula, standard_material, assemblage):
    """
        Parameters
        ----------
        component_formula : dictionary
            Chemical formula dictionary

        standard_material: class
            Material class
            set_method and set_state should already have been used

        assemblage: list of classes
            List of material classes
            set_method and set_state should already have been used

        Returns
        -------
        fugacity : float
            Value of the fugacity of the component with respect to
            the standard material

    """
    chemical_potential=chemicalpotentials(assemblage, [component_formula])[0]

    fugacity=np.exp((chemical_potential - standard_material.gibbs)/(R*assemblage[0].temperature))
    return fugacity

def relativefugacity(component_formula, standard_material, assemblage, reference_assemblage):
    """
        Parameters
        ----------
        component_formula : dictionary
            Chemical formula dictionary

        standard_material: class
            Material class
            set_method and set_state should already have been used

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
    chemical_potential=chemicalpotentials(assemblage, [component_formula])[0]
    reference_chemical_potential=chemicalpotentials(reference_assemblage, [component_formula])[0]

    relative_fugacity=np.exp((chemical_potential - reference_chemical_potential)/(R*assemblage[0].temperature))
    return relative_fugacity
