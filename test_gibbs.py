import burnman
from burnman.minerals.SLB_2011 import *
import numpy as np
import matplotlib.pyplot as plt


minlist = [mg_fe_olivine(), wuestite(), periclase(), stishovite(), jadeite()]
composition = {'Mg': 0.2, 'Si': 0.2, 'Fe': 0.1, 'O': 0.5}

ea = burnman.EquilibriumAssemblage(composition, minlist)
ea.set_method('slb3')
ea.set_state(100.e9, 3000.)
