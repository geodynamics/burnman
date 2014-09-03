import os, sys, numpy as np, matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
from chemicalpotentials import *
import numpy as np

'''
Here we initialise the minerals we'll be using
'''
P=1.e8
T=1000.

py=minerals.HP_2011_ds62.py()
fo=minerals.HP_2011_ds62.fo()
en=minerals.HP_2011_ds62.en()
per=minerals.HP_2011_ds62.per()
minerals=[py, fo, en, per]

for mineral in minerals:
    mineral.set_method('mtait')
    mineral.set_state(P, T)


'''
Here we find chemical potentials of MgO, Al2O3 and SiO2 for
an assemblage containing pyrope, forsterite and enstatite 
at 0.1 GPa, 1000 K
'''

assemblage=[py, fo, en]
component_formulae=['MgO', 'Al2O3', 'SiO2']
component_formulae_dict=[dictionarize_formula(f) for f in component_formulae]
chem_potentials=chemicalpotentials(assemblage, component_formulae_dict)

for idx, component_formula in enumerate(component_formulae):
    print "mu("+component_formula+"):", chem_potentials[idx]
print ''


'''
Here we find the 'periclase fugacity' of the assemblage
Fugacity is often defined relative to some reference pressure
Here we use room pressure, 100 kPa
'''

Pr=1.e5
per.set_state(Pr, T)
component=dictionarize_formula('MgO')
print 'log(fugacity_per):', np.log10(fugacity(component, per, assemblage))
