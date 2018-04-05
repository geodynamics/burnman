# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.optimize import fsolve
import burnman

liq = burnman.minerals.other.liquid_iron()
burnman.tools.check_eos_consistency(liq, P=10.e9, T=7000., tol=1.e-3, verbose=True)


# Find heat capacities
temperatures = np.linspace(1000., 15000., 101)
Cvs = np.empty_like(temperatures)
m = 0.055845
rhos = np.empty_like(temperatures)
densities = [5.e3,10.e3, 15.e3]
for rho in densities:
    V = m/rho
    for i, T in enumerate(temperatures):
        Cvs[i] = liq.method.molar_heat_capacity_v(0., T, V, liq.params)/burnman.constants.gas_constant

    plt.plot(temperatures, Cvs)

fig1 = mpimg.imread('../../burnman/data/input_figures/AA1994_liq_iron_TCv_different_densities.png')
plt.imshow(fig1, extent=[1000., 15000., 0., 6.], aspect='auto')
plt.ylim(0., 6.)
plt.title('AA1994, Figure 5')
plt.xlabel('Temperature (K)')
plt.ylabel('Cv/R/atom')
plt.show()


# Find volumes and temperatures up the reference isentrope
# Check standard state values
liq.set_state(1.e5, 1811.)
reference_entropy = liq.S

def isentrope(T, P, Sref, mineral):
    mineral.set_state(P, T[0])
    return Sref - mineral.S

pressures = np.linspace(0., 500.e9, 21)
temperatures = np.empty_like(pressures)
rhos = np.empty_like(pressures)
Vps = np.empty_like(pressures)

for i, P in enumerate(pressures):
    temperatures[i] = fsolve(isentrope, 1811., args=(P, reference_entropy, liq))[0]
    rhos[i] = liq.density

fig1 = mpimg.imread('../../burnman/data/input_figures/AA1994_liq_iron_PTrho_reference_isentrope.png')

plt.imshow(fig1, extent=[0.0, 500., 6., 15.], aspect='auto')
plt.plot(pressures/1.e9, rhos/1.e3, marker='o', linestyle='None') 
plt.title('1811 K isentrope; AA1994 Figure B1 (1/2)')
plt.xlabel('Pressure (GPa)')
plt.ylabel('Density (kg/m^3)')
plt.show()

plt.imshow(fig1, extent=[0.0, 500., 1500., 7000.], aspect='auto')
plt.plot(pressures/1.e9, temperatures, marker='o', linestyle='None')
plt.title('1811 K isentrope; AA1994 Figure B1 (2/2)')
plt.xlabel('Pressure (GPa)')
plt.ylabel('Temperature (K)')
plt.show()


# Find densities at 1 bar
temperatures = np.linspace(1800., 2400., 100)
rhos = np.empty_like(temperatures)
rhos_mizuno = np.empty_like(temperatures)

P = 1.e5
for i, T in enumerate(temperatures):
    liq.set_state(1.e5, T)
    rhos[i] = liq.density/1.e3
    rhos_mizuno[i] = (7162. - 0.735*(T - 1808.))/1.e3



fig1 = mpimg.imread('../../burnman/data/input_figures/AA1994_liq_iron_Trho_1bar.png')
plt.imshow(fig1, extent=[1800., 2400., 6.65, 7.1], aspect='auto')
plt.plot(temperatures, rhos, label='Model')
plt.plot(temperatures, rhos_mizuno, label='Mizuno et al.')
plt.ylim(6.65, 7.1)
plt.xlabel('Temperatures (K)')
plt.ylabel("Density (kg/m^3)")
plt.title('1 bar densities; AA1994 Figure 1')
plt.legend(loc='lower left')
plt.show()


def temperature(T, P, rho, mineral):
    mineral.set_state(P, T[0])
    return mineral.density - rho

# Check grueneisen values
Prhos = [[1.e5, 7019.],
         [0.2e9, 5500.],
         [0.2e9, 6000.],
         [0.2e9, 6500.],
         [277.4e9, 12643.],
         [331.5e9, 13015.],
         [397.1e9, 13417.]]

print('Pressure (GPa), Temperature (K), Density (kg/m^3), Grueneisen parameter')
for Prho in Prhos:
    P, rho = Prho
    T = fsolve(temperature, 1811., args=(P, rho, liq))[0]
    print('{0:.1f} {1:.0f} {2:.1f} {3:.4f}'.format(P/1.e9, rho, T, liq.grueneisen_parameter))
