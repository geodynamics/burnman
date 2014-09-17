import burnman
from burnman.minerals.SLB_2011 import *
import numpy as np
import matplotlib.pyplot as plt


minlist = [mg_fe_olivine(), wuestite(), periclase(), stishovite()]

M, elements, formulae = burnman.assemble_stoichiometric_matrix( minlist )
print M
N = burnman.compute_nullspace(M)
print N
