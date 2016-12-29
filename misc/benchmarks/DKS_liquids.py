import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman.minerals import DKS_2013_liquids
from burnman import constants
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


SiO2_liq=DKS_2013_liquids.SiO2_liquid()
MgO_liq=DKS_2013_liquids.MgO_liquid()

fig1 = mpimg.imread('figures/SiO2_liquid_PVT.png')
plt.imshow(fig1, extent=[9, 30, -10, 220], aspect='auto')

temperatures=np.linspace(2000., 7000., 6)
volumes=np.linspace(9e-6, 30e-6, 101)
pressures=np.empty_like(volumes)
for temperature in temperatures:
    for i, volume in enumerate(volumes):
        pressures[i]=SiO2_liq.method.pressure(temperature, volume, SiO2_liq.params)/1e9
    plt.plot(volumes*1e6, pressures, linewidth=2, label=str(temperature)+'K')

plt.legend(loc='upper right')
plt.xlim(9,30)
plt.ylim(-10,220)
plt.show()

fig1 = mpimg.imread('figures/SiO2_liquid_SelVT.png')
plt.imshow(fig1, extent=[9, 30, -0.03, 0.75], aspect='auto')

temperatures=np.linspace(2000., 7000., 6)
volumes=np.linspace(9e-6, 30e-6, 101)
entropies=np.empty_like(volumes)
for temperature in temperatures:
    for i, volume in enumerate(volumes):
        entropies[i]=SiO2_liq.method._S_el(temperature, volume, SiO2_liq.params)
    plt.plot(volumes*1e6, entropies/3./constants.gas_constant, linewidth=2, label=str(temperature)+'K')

plt.legend(loc='upper right')
plt.xlim(9,30)
plt.show()

fig1 = mpimg.imread('figures/SiO2_liquid_EVT.png')
plt.imshow(fig1, extent=[9, 30, -2400, -1200], aspect='auto')

temperatures=np.linspace(2000., 7000., 6)
volumes=np.linspace(9e-6, 30e-6, 101)
energies=np.empty_like(volumes)
for temperature in temperatures:
    for i, volume in enumerate(volumes):
        energies[i]=SiO2_liq.method.internal_energy(0., temperature, volume, SiO2_liq.params)
    plt.plot(volumes*1e6, energies/1e3, linewidth=2, label=str(temperature)+'K')

plt.legend(loc='upper right')
plt.xlim(9,30)
plt.show()

# Plot MgO
fig1 = mpimg.imread('figures/MgO_liquid_PVT.png')
plt.imshow(fig1, extent=[7, 18, -6, 240], aspect='auto')

temperatures=np.linspace(2000., 7000., 6)
volumes=np.linspace(7e-6, 18e-6, 101)
pressures=np.empty_like(volumes)
for temperature in temperatures:
    for i, volume in enumerate(volumes):
        pressures[i]=SiO2_liq.method.pressure(temperature, volume, MgO_liq.params)/1e9
    plt.plot(volumes*1e6, pressures, linewidth=2, label=str(temperature)+'K')

temperature = 10000. # K
for i, volume in enumerate(volumes):
    pressures[i]=SiO2_liq.method.pressure(temperature, volume, MgO_liq.params)/1e9
plt.plot(volumes*1e6, pressures, linewidth=2, label=str(temperature)+'K')

    
plt.legend(loc='upper right')
plt.xlim(7,18)
plt.ylim(-6,240)
plt.show()


fig1 = mpimg.imread('figures/MgO_liquid_SelVT.png')
plt.imshow(fig1, extent=[6, 18, -0.04, 0.84], aspect='auto')

temperatures=np.linspace(2000., 7000., 6)
volumes=np.linspace(7e-6, 18e-6, 101)
entropies=np.empty_like(volumes)
for temperature in temperatures:
    for i, volume in enumerate(volumes):
        entropies[i]=SiO2_liq.method._S_el(temperature, volume, MgO_liq.params)
    plt.plot(volumes*1e6, entropies/2./constants.gas_constant, linewidth=2, label=str(temperature)+'K')

temperature = 10000. # K
for i, volume in enumerate(volumes):
    entropies[i]=SiO2_liq.method._S_el(temperature, volume, MgO_liq.params)
plt.plot(volumes*1e6, entropies/2./constants.gas_constant, linewidth=2, label=str(temperature)+'K')

    
plt.legend(loc='upper right')
plt.xlim(6,18)
plt.ylim(-0.04,0.84)
plt.show()


fig1 = mpimg.imread('figures/MgO_liquid_EVT.png')
plt.imshow(fig1, extent=[7, 18, -1200, -200], aspect='auto')

temperatures=np.linspace(2000., 7000., 6)
volumes=np.linspace(7e-6, 18e-6, 101)
energies=np.empty_like(volumes)
for temperature in temperatures:
    for i, volume in enumerate(volumes):
        energies[i]=MgO_liq.method.internal_energy(0., temperature, volume, MgO_liq.params)
    plt.plot(volumes*1e6, energies/1e3, linewidth=2, label=str(temperature)+'K')

temperature=10000. # K
for i, volume in enumerate(volumes):
    energies[i]=MgO_liq.method.internal_energy(0., temperature, volume, MgO_liq.params)
plt.plot(volumes*1e6, energies/1e3, linewidth=2, label=str(temperature)+'K')

plt.legend(loc='upper right')
plt.xlim(7,18)
plt.show()
