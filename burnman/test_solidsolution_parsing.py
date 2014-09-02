# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# This is a test script which is the precursor to implementing a solid solution BaseClass.


'''
# Inputs
Solid solution model
Endmember proportions
P,T

# Things we can initialise without endmember proportions
n_sites
n_occupancies
n_endmembers
alpha
Wh, Ws, Wv
site_occupancies
site_multiplicities

# Endmember proportions needed
phi
ideal activities
occupancies

# P, T needed
Wtotal / nonideal contribution to gibbs
'''

import re
import numpy as np
from fractions import Fraction
from processchemistry import ProcessSolidSolutionChemistry

# Endmembers
base_material = [['pyrope()',    '[Mg]3[Al]2Si3O12',         1.0], \
                 ['almandine()', '[Fe]3[Al]2Si3O12',         1.0], \
                 ['grossular()', '[Ca]3[Al]2Si3O12',         2.7], \
                 ['majorite()',  '[Mg]3[Mg1/2Si1/2]2Si3O12', 1.0]]

n_endmembers=len(base_material)

# Interaction parameters
excess_enthalpy=[[2.5e3, 30.1e3, 15e3],[10e3,18e3],[48e3]]
excess_entropy=[[0., 0., 0.],[0., 0.],[0.]]
excess_volume=[[0., 0.164e-5, 0.],[0., 0.],[0.]]
interaction_parameter=[excess_enthalpy,excess_entropy,excess_volume]


# INPUT PROPORTIONS
molar_fraction = np.array([ 0.5, 0.2, 0.1, 0.2 ])



solution_formulae, n_sites, sites, n_occupancies, endmember_occupancies, site_multiplicities = ProcessSolidSolutionChemistry([base_material[i][1] for i in range(len(base_material))])


print solution_formulae
print n_sites
print n_occupancies
print sites
print endmember_occupancies
print site_multiplicities 

