import burnman
from burnman.minerals.SLB_2011 import *
import numpy as np
import matplotlib.pyplot as plt


composition = { 'Mg': 2./7., 'O': 4./7., 'Si': 1./7}
minlist = [forsterite(), mg_wadsleyite(), mg_ringwoodite(), mg_perovskite(), periclase(), stishovite(),mg_akimotoite(), mg_majorite() ]

fo = forsterite()
fo.set_method('slb3')
wa = mg_wadsleyite()
wa.set_method('slb3')
ri = mg_ringwoodite()
ri.set_method('slb3')
pvpe = burnman.Composite( [0.5, 0.5], [mg_perovskite(), periclase()])
pvpe.set_method('slb3')
pest = burnman.Composite( [0.5, 0.5], [periclase(), stishovite()])
pest.set_method('slb3')

ea = burnman.EquilibriumAssemblage(composition, minlist)
ea.set_method('slb3')

n= 100.
pressures = np.linspace(1.e5, 30.e9, n)
temperatures = np.linspace(800., 2800.0, n)
G = np.empty_like(pressures)
Gr = np.empty_like(pressures)
Gf = np.empty_like(pressures)
Gw = np.empty_like(pressures)
Gp = np.empty_like(pressures)
Gs = np.empty_like(pressures)
Ga = np.empty_like(pressures)

phase = np.empty( (n,n) )
T = 1600.

for i,p in enumerate(pressures):
    ea.set_state(p, T)
    G[i] = ea.gibbs*7.
    fo.set_state(p, T)
    Gf[i] = fo.gibbs
    wa.set_state(p, T)
    Gw[i] = wa.gibbs
    ri.set_state(p, T)
    Gr[i] = ri.gibbs
    pvpe.set_state(p, T)
    Gp[i] = pvpe.gibbs*2.
   

Ga = (Gw+Gr+Gf+Gp)*1./4.
plt.plot(pressures, G-Ga, 'k--', linewidth=4)
plt.plot(pressures, Gr-Ga)
plt.plot(pressures, Gw-Ga)
plt.plot(pressures, Gf-Ga)
plt.plot(pressures, Gp-Ga)
plt.show()
plt.clf()

eps = 1.e-4
for i,p in enumerate(pressures):
    for j,t in enumerate(temperatures):
        ea.set_state(p, t)
        #forsterite field
        if ea.species_vector[0]/np.sum(ea.species_vector) > eps:
            phase[j,i] = 0.
        #wadsleyite field
        if ea.species_vector[1]/np.sum(ea.species_vector) > eps:
            phase[j,i] = 1.
        #ringwoodite field
        if ea.species_vector[2]/np.sum(ea.species_vector) > eps:
            phase[j,i] = 2.
        #perovskite+periclase field
        if ea.species_vector[3]/np.sum(ea.species_vector) > eps and ea.species_vector[4]/np.sum(ea.species_vector) > eps:
            phase[j,i] = 3.
        #periclase + stishovite field
        if ea.species_vector[4]/np.sum(ea.species_vector) > eps and ea.species_vector[5]/np.sum(ea.species_vector) > eps:
            phase[j,i] = 4.
        #akimotoite + periclase field
        if ea.species_vector[4]/np.sum(ea.species_vector) > eps and ea.species_vector[6]/np.sum(ea.species_vector) > eps:
            phase[j,i] = 5.
        #majorite + periclase field
        if ea.species_vector[4]/np.sum(ea.species_vector) > eps and ea.species_vector[7]/np.sum(ea.species_vector) > eps:
            phase[j,i] = 6.

plt.imshow( phase, origin='lower', extent=[0., 30., 800., 2800.], aspect='auto')
plt.xlabel("Pressure (GPa)")
plt.ylabel("Temperature (K)")
plt.show()
