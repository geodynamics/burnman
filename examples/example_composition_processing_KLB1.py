from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))


import burnman
from burnman.processchemistry import dictionarize_formula

# Let's import a function to read compositions from file, and the Composition class
# The function assumes that the first line is a header of components, and then the word "Comment"
# The following lines are either molar or weight proportions of the components given in the header,
# followed by any string
from burnman.composition import Composition

# We want to remove a fraction of fo90 from our starting composition
# Here, I chose an amount that would produce a bulk composition that was a lot less
# olivine-polymorph rich, but still contained a "safe" amount of the polymorph (ca. 20%)
# (for buffering purposes).
# The exact amount is not important, 35.8% (oxide basis) produced exactly 20% ol given the
# composition of the majorite from one of the 14 GPa runs.
x = 35.8 
KLB1 = Composition({'SiO2': 39.4 - x*1./3.,
                    'Al2O3': 2.0,
                    'CaO': 3.3,
                    'MgO': 49.5 - x*2./3.*0.9,
                    'FeO': 5.2 - x*2./3.*0.1}, 'molar') # from Holland et al., 2013



# Now we need to modify the composition so that we can make it from
# standard starting materials:

# 1) We want to add calcium via CaCO3.

# 2) We want to make a fraction f of total Fe Fe57 (where 1 is 100%)
# where the Fe57 is metallic Fe and the Fe(natural) is Fe2O3
# Obviously we only want to vary the amount of O2, so
# FeO + m O -> f Fe(57) + (1-f)/2 Fe2O3

# Therefore m = 3*(1-f)/2 - 1 = 0.5 - 1.5f
f = 0.
m = 0.5 - 1.5*f

# Here's where we add the (temporary) components 
KLB1.add_components(composition_dictionary = {'CO2': KLB1.molar_composition['CaO'],
                                              'O': KLB1.molar_composition['FeO']*m},
                    unit_type='molar')

# Here's where we change the components:
KLB1.change_component_set(['CaCO3', 'Fe', 'Fe2O3', 'MgO', 'Al2O3', 'SiO2'])
KLB1.print('weight', significant_figures=4, total=1.)


