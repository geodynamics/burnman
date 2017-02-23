from __future__ import absolute_import
from __future__ import print_function
import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman.minerals import \
    DKS_2013_liquids, \
    DKS_2013_solids, \
    SLB_2011
from burnman import constants
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

phases = [DKS_2013_liquids.SiO2_liquid(),
         DKS_2013_liquids.MgSiO3_liquid(),
         DKS_2013_liquids.MgSi2O5_liquid(),
         DKS_2013_liquids.MgSi3O7_liquid(),
         DKS_2013_liquids.MgSi5O11_liquid(),
         DKS_2013_liquids.Mg2SiO4_liquid(),
         DKS_2013_liquids.Mg3Si2O7_liquid(),
         DKS_2013_liquids.Mg5SiO7_liquid(),
         DKS_2013_liquids.MgO_liquid()
         ]


pressure = 25.e9 # Pa
temperature = 3000. # K

MgO_liq = DKS_2013_liquids.MgO_liquid()
SiO2_liq = DKS_2013_liquids.SiO2_liquid()

pressures = np.array(   [0.,    5.,    10.,   25.,   50.,   75.,   100.,  135.])*1.e9
temperatures = np.array([3000., 3000., 3000., 3000., 4000., 5000., 6000., 6000.])

SiO2_enthalpies, SiO2_entropies, SiO2_volumes = phases[0].evaluate(['H', 'S', 'V'], pressures, temperatures)
MgO_enthalpies, MgO_entropies, MgO_volumes = phases[-1].evaluate(['H', 'S', 'V'], pressures, temperatures)


excesses = []
for phase in phases[1:-1]:
    try:
        nSi = phase.params['formula']['Si']
    except:
        nSi = 0.
    try:
        nMg = phase.params['formula']['Mg']
    except:
        nMg = 0.
        
    sum_cations = nSi+nMg
    fSi=nSi/sum_cations

    enthalpies, entropies, volumes = phase.evaluate(['H', 'S', 'V'], pressures, temperatures)

    excess_enthalpies = enthalpies/sum_cations - (fSi*SiO2_enthalpies + (1.-fSi)*MgO_enthalpies)
    excess_entropies = entropies/sum_cations - (fSi*SiO2_entropies + (1.-fSi)*MgO_entropies)
    excess_volumes = volumes/sum_cations - (fSi*SiO2_volumes + (1.-fSi)*MgO_volumes)



    excesses.extend([[pressure, temperatures[i], fSi,
                      excess_enthalpies[i], excess_entropies[i], excess_volumes[i]]
                     for i, pressure in enumerate(pressures)])


excesses = np.array(excesses)

fig = plt.figure()
ax_H = fig.add_subplot(1,3,1)
ax_S = fig.add_subplot(1,3,2)
ax_V = fig.add_subplot(1,3,3)

ax_H.set_title('Excess Enthalpy') 
ax_S.set_title('Excess Entropy')
ax_V.set_title('Excess Volume') 

figH = mpimg.imread('figures/MgO-SiO2_enthalpy_mixing.png')
figS = mpimg.imread('figures/MgO-SiO2_entropy_mixing.png')
figV = mpimg.imread('figures/MgO-SiO2_volume_mixing.png')

ax_H.imshow(figH, extent=[0, 1, -40000, 10000], aspect='auto')
ax_S.imshow(figS, extent=[0, 1, 0, 24], aspect='auto')
ax_V.imshow(figV, extent=[0, 1, -2.0e-6, 0.4e-6], aspect='auto')

for (pressure, temperature) in zip(*[pressures, temperatures]):
    mask = [i for i in range(len(excesses[:,1])) if np.abs(excesses[i][0] - pressure) < 1. and np.abs(excesses[i][1] - temperature) < 1.]
    ax_H.plot(excesses[mask,2], excesses[mask,3], marker='o', linestyle='None', label='{0:.0f} GPa, {1:.0f} K'.format(pressure/1.e9, temperature))
    ax_S.plot(excesses[mask,2], excesses[mask,4], marker='o', linestyle='None', label='{0:.0f} GPa, {1:.0f} K'.format(pressure/1.e9, temperature))
    ax_V.plot(excesses[mask,2], excesses[mask,5], marker='o', linestyle='None', label='{0:.0f} GPa, {1:.0f} K'.format(pressure/1.e9, temperature))

plt.legend(loc='lower right')
plt.show()
    
print('Excess properties of mixing at a 1:1 Mg:Si composition:\n')
mask = [i for i in range(len(excesses[:,2])) if np.abs(excesses[i][2] - 0.5) < 1.e-12]
print('P (GPa)   T (K)   H (J/mol)   S (J/K/mol)   V (cm^3/mol)')
for excess in excesses[mask]:
    print('{0:.0f} {1:.0f} {2:.1e} {3:.1e} {4:.1e}'.format(excess[0]/1.e9, excess[1], excess[3], excess[4], excess[5]))


