# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import burnman
from burnman import minerals
from processchemistry import *
import numpy as np

# This module computes chemical potentials (partial molar gibbs free energies) for an assemblage based on the Gibbs free energies and compositions of the individual phases.

# It can also calculate fugacities based on the gibbs free energies of the endmembers corresponding to chemical components.

# It can also calculate fugacities relative to other bulk compositions

def chemicalpotentials(minerals, component_formulae):
    mineral_gibbs=[]
    mineral_formulae=[]
    for mineral in minerals:
        mineral_gibbs.append(mineral.gibbs)
        mineral_formulae.append(mineral.params['formula'])
        
    mineral_compositions, elements=compositional_array(mineral_formulae)

    for i in range(len(component_formulae)):
        component_formulae[i]=dictionarize_formula(component_formulae[i])

    component_compositions=ordered_compositional_array(component_formulae, elements)
    component_proportions=np.linalg.lstsq(component_compositions.T, mineral_compositions.T)[0].T
    component_potentials=np.linalg.solve(component_proportions,mineral_gibbs)

    return component_potentials

