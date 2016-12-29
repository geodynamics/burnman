import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman.minerals import DKS_2013_solids, DKS_2008_fo
from burnman import constants
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

'''
SOLIDS
'''

# Vxfo = fo.params['V_x']*1.e6
#['forsterite', DKS_2008_fo.forsterite(), [0.67*Vxfo, 1.03*Vxfo, -15, 55], [0.67*Vxfo, 1.03*Vxfo, -2400, -1800]], \
#['fo_liquid', DKS_2008_fo.fo_liquid(), [0.48*Vxfo, 1.23*Vxfo, -5, 205], [0.48*Vxfo, 1.23*Vxfo, -2400, -1800]], \
          
phases = [['stishovite', DKS_2013_solids.stishovite(), [10, 18, -25, 175], [10, 18, -2400, -1800]], \
          ['perovskite', DKS_2013_solids.perovskite(), [14.5, 27.5, 0, 344], [14.5, 27.5, -3600, -2000]], \
          ['periclase', DKS_2013_solids.periclase(), [6.5, 14, -25, 275], [6.5, 14, -1200, -560]]]

for name, phase, PVT_range, EVT_range in phases:
    print name
    vmin=PVT_range[0]
    vmax=PVT_range[1]
 
    fig1 = mpimg.imread('figures/'+name+'_PVT.png')
    plt.imshow(fig1, extent=PVT_range, aspect='auto')
    
    temperatures=np.linspace(1000., 6000., 6)
    volumes=np.linspace(PVT_range[0]*1.e-6, PVT_range[1]*1.e-6, 101)
    pressures=np.empty_like(volumes)
    for temperature in temperatures:
        for i, volume in enumerate(volumes):
            pressures[i]=phase.method.pressure(temperature, volume, phase.params)
        plt.plot(volumes*1e6, pressures/1e9, linewidth=2, label=str(temperature)+'K')

    plt.legend(loc='upper right')
    plt.xlim(PVT_range[0], PVT_range[1])
    plt.show()
    


    fig1 = mpimg.imread('figures/'+name+'_EVT.png')
    plt.imshow(fig1, extent=EVT_range, aspect='auto')
    
    temperatures=np.linspace(2000., 6000., 5)
    volumes=np.linspace(EVT_range[0]*1.e-6, EVT_range[1]*1.e-6, 101)
    energies=np.empty_like(volumes)
    for temperature in temperatures:
        for i, volume in enumerate(volumes):
            energies[i]=phase.method.internal_energy(0., temperature, volume, phase.params)
        plt.plot(volumes*1e6, energies/1e3, linewidth=2, label=str(temperature)+'K')

    plt.legend(loc='upper right')
    plt.xlim(EVT_range[0], EVT_range[1])
    plt.show()
