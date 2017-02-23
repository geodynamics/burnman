from __future__ import absolute_import
from __future__ import print_function
import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman.minerals import DKS_2013_liquids
from burnman import constants
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


SiO2_liq=DKS_2013_liquids.SiO2_liquid()
MgO_liq=DKS_2013_liquids.MgO_liquid()

MgO_temperatures=np.array([2000., 3000., 4000., 5000., 6000., 7000., 10000.])
MgO_volumes=np.linspace(7e-6, 18e-6, 101)

SiO2_temperatures=np.array([2000., 3000., 4000., 5000., 6000., 7000.])
SiO2_volumes=np.linspace(9e-6, 30e-6, 101)

plots = [[['figures/SiO2_liquid_PVT.png', [9, 30, -10, 220]],
          ['figures/SiO2_liquid_SelVT.png', [9, 30, -0.03, 0.75]],
          ['figures/SiO2_liquid_EVT.png', [9, 30, -2400, -1200]]],
         [['figures/MgO_liquid_PVT.png', [7, 18, -6, 240]],
          ['figures/MgO_liquid_SelVT.png', [6, 18, -0.04, 0.84]],
          ['figures/MgO_liquid_EVT.png', [7, 18, -1200, -200]]]]


for i, (phase, n_atoms, temperatures, volumes) in enumerate([(SiO2_liq, 3., SiO2_temperatures,
                                                              SiO2_volumes),
                                                             (MgO_liq, 2., MgO_temperatures,
                                                              MgO_volumes)]):

    fig = plt.figure()
    ax_P = fig.add_subplot(1,3,1)
    ax_S = fig.add_subplot(1,3,2)
    ax_E = fig.add_subplot(1,3,3)
    ax_P.set_xlabel('Volume (cm^3/mol)')
    ax_S.set_xlabel('Volume (cm^3/mol)')
    ax_E.set_xlabel('Volume (cm^3/mol)')
    
    ax_P.set_ylabel('Pressure (GPa)')
    ax_S.set_ylabel('Entropy/nR (J)')
    ax_E.set_ylabel('Internal Energy (kJ/mol)')

    ax_P.imshow(mpimg.imread(plots[i][0][0]), extent=plots[i][0][1], aspect='auto')
    ax_S.imshow(mpimg.imread(plots[i][1][0]), extent=plots[i][1][1], aspect='auto')
    ax_E.imshow(mpimg.imread(plots[i][2][0]), extent=plots[i][2][1], aspect='auto')
    
    for temperature in temperatures:
        pressures = np.empty_like(volumes)
        entropies= np.empty_like(volumes)
        energies = np.empty_like(volumes)

        for j, volume in enumerate(volumes):
            pressures[j]=phase.method.pressure(temperature, volume, phase.params)
            entropies[j]=phase.method._S_el(temperature, volume, phase.params)
            energies[j]=phase.method.internal_energy(0., temperature, volume, phase.params)

        ax_P.plot(volumes*1.e6, pressures/1.e9, linewidth=2, label='{0:.0f} K'.format(temperature))
        ax_S.plot(volumes*1.e6, entropies/n_atoms/constants.gas_constant, linewidth=2, label='{0:.0f} K'.format(temperature))
        ax_E.plot(volumes*1.e6, energies/1.e3, linewidth=2, label='{0:.0f} K'.format(temperature))
        
        ax_E.legend(loc='upper right')
    plt.show()
