import burnman
from burnman.minerals.SLB_2011 import *
import numpy as np
import matplotlib.pyplot as plt


minlist = [mg_fe_olivine(), ferropericlase(), stishovite()]
composition = {'Fe': 1./7., 'Mg': 1./7., 'O': 4./7., 'Si': 1./7}

ea = burnman.EquilibriumAssemblage(composition, minlist)
ea.set_method('slb3')
