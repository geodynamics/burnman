import burnman
from burnman.minerals.SLB_2011 import *
import numpy as np
import matplotlib.pyplot as plt


#minlist = [mg_fe_olivine(), mg_fe_wadsleyite()]
#composition = {'Fe': 1./35., 'Mg': 9./35., 'O': 4./7., 'Si': 1./7}
composition = { 'Mg': 2./7., 'O': 4./7., 'Si': 1./7}
minlist = [forsterite(), mg_wadsleyite(), mg_ringwoodite(), mg_perovskite(), periclase(), stishovite()]
#minlist = [forsterite(), mg_wadsleyite(), mg_ringwoodite()]

fo = forsterite()
fo.set_method('slb3')
wa = mg_wadsleyite()
wa.set_method('slb3')



ea = burnman.EquilibriumAssemblage(composition, minlist)
ea.set_method('slb3')

n= 60.
pressures = np.linspace(1.e5, 30.e9, n)
temperatures = np.linspace(800., 2800.0, n)

phase = np.empty( (n,n) )
eps = 1.e-5



#for i,p in enumerate(pressures):
#    for j,t in enumerate(temperatures):
#        ea.set_state(p, t)
#        phase[j,i] = np.argmax(ea.species_vector)
#plt.imshow( np.flipud(phase))
#plt.colorbar()
#plt.show()
        
