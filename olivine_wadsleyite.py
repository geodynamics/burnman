import burnman
from burnman.minerals.SLB_2011 import *
import numpy as np
import matplotlib.pyplot as plt


#composition = { 'Fe': 1./35., 'Mg': 9./35., 'O': 4./7., 'Si': 1./7.}
composition = { 'Mg': 1.8, 'Fe': 0.2, 'O': 4., 'Si': 1.}
minlist = [mg_fe_olivine(), mg_fe_wadsleyite()]

ol = mg_fe_olivine()
ol.set_method('slb3')
ol.set_composition( [0.9, 0.1])
wa = mg_fe_wadsleyite()
wa.set_method('slb3')
wa.set_composition( [0.9, 0.1])
ea = burnman.EquilibriumAssemblage(composition, minlist)
ea.set_method('slb3')
co = burnman.Composite( [0.5, 0.5], [ol, wa] )

pressures = np.linspace( 12.e9, 14.e9, 50 )
temperature = 1500.
Gol = np.empty_like(pressures)
Gwa = np.empty_like(pressures)
Gea = np.empty_like(pressures)
Gav = np.empty_like(pressures)


for i,p in enumerate(pressures):
    ol.set_state(p,temperature)
    wa.set_state(p,temperature)
    ea.set_state(p,temperature)
    ea.print_assemblage()

    Gol[i] = ol.gibbs
    Gwa[i] = wa.gibbs
    Gea[i] = ea.gibbs


Gav = (Gol+Gwa)*1./2.
plt.plot(pressures, Gol-Gav)
plt.plot(pressures, Gwa-Gav)
plt.plot(pressures, Gea-Gav, 'k-', linewidth=2)
plt.show()
