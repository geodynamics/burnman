import os, sys, numpy as np, matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
from chemicalpotentials import *
import numpy as np


P=1.e8
T=1000.

py=minerals.HP_2011_ds62.py()
fo=minerals.HP_2011_ds62.fo()
en=minerals.HP_2011_ds62.en()

minerals=[py, fo, en]

for mineral in minerals:
    mineral.set_method('mtait')
    mineral.set_state(P, T)

component_formulae=['MgO', 'Al2O3', 'SiO2']

chem_potentials=chemicalpotentials(minerals, component_formulae)

for idx, mineral in enumerate(minerals):
    print mineral.params['name'], chem_potentials[idx]

