from __future__ import absolute_import
from __future__ import print_function
import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman.minerals import DKS_2013_solids
from burnman import constants
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

'''
SOLIDS
'''
          
phases = [['stishovite', DKS_2013_solids.stishovite(), [10, 18, -25, 175], [10, 18, -2400, -1800]], \
          ['perovskite', DKS_2013_solids.perovskite(), [14.5, 27.5, 0, 344], [14.5, 27.5, -3600, -2000]], \
          ['periclase', DKS_2013_solids.periclase(), [6.5, 14, -25, 275], [6.5, 14, -1200, -560]]]

temperatures=[1000., 2000., 3000., 4000., 5000., 6000., 8000.]
for name, phase, PVT_range, EVT_range in phases:

    phase.params['G_0'] = 0. # just for consistency checking
    phase.params['Gprime_0'] = 1.3 # just for consistency checking
    phase.params['eta_s_0'] = 0. # just for consistency checking
    print('EoS consistent for {0} model: {1}'.format(name, burnman.tools.check_eos_consistency(phase, tol=1.e-4)))
    
    vmin=PVT_range[0]
    vmax=PVT_range[1]

    fig = plt.figure()
    ax_P = fig.add_subplot(1, 2, 1)
    ax_E = fig.add_subplot(1, 2, 2)
    
    fig1 = mpimg.imread('figures/'+name+'_PVT.png')
    ax_P.imshow(fig1, extent=PVT_range, aspect='auto')
    
    volumes=np.linspace(PVT_range[0]*1.e-6, PVT_range[1]*1.e-6, 101)
    pressures=np.empty_like(volumes)
    for temperature in temperatures:
        for i, volume in enumerate(volumes):
            pressures[i]=phase.method.pressure(temperature, volume, phase.params)
        ax_P.plot(volumes*1e6, pressures/1e9, linewidth=2, label='{0:.0f} K'.format(temperature))

    ax_P.set_xlim(PVT_range[0], PVT_range[1])
    ax_P.set_xlabel('Volume (cm^3/mol)')
    ax_P.set_ylabel('Pressure (GPa)')
    

    fig1 = mpimg.imread('figures/'+name+'_EVT.png')
    ax_E.imshow(fig1, extent=EVT_range, aspect='auto')

    volumes=np.linspace(EVT_range[0]*1.e-6, EVT_range[1]*1.e-6, 101)
    energies=np.empty_like(volumes)
    for temperature in temperatures:
        for i, volume in enumerate(volumes):
            energies[i]=phase.method.molar_internal_energy(0., temperature, volume, phase.params)
        ax_E.plot(volumes*1e6, energies/1e3, linewidth=2, label='{0:.0f} K'.format(temperature))

    ax_E.legend(loc='upper right')
    ax_E.set_xlim(EVT_range[0], EVT_range[1])
    ax_E.set_xlabel('Volume (cm^3/mol)')
    ax_E.set_ylabel('Internal energy (kJ/mol)')
    
    fig.canvas.set_window_title(name)
    plt.show()
