# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

example_thermodynamic_beginner
----------------------
    
This example is for absolute beginners wishing to use burnman to do 
thermodynamics. 

*Uses:*

* :doc:`mineral_database`


*Demonstrates:*

* Ways to query endmember properties
* How to set a reference temperature

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
import matplotlib.pyplot as plt

fo=minerals.HP_2011_ds62.fo()
fo.params['P_0'] = 0.99999999999e5
room_pressure = 1.e5
temperatures = np.linspace(50., 1500., 100)
alphas_0 = []
volumes_0 = []
for T in temperatures:
    fo.set_state(room_pressure, T)
    alphas_0.append(fo.alpha)
    volumes_0.append(fo.V)

T_eos = 2273.
pressures = np.linspace(1.e5, 50.e9, 100)
eos_0 = []
for P in pressures:
    fo.set_state(P, T_eos)
    eos_0.append(fo.V)


T_ref = 1273.
fo.set_state(room_pressure, T_ref)
fo.params['T_0'] = T_ref
fo.params['H_0'] = fo.H
fo.params['S_0'] = fo.S
fo.params['V_0'] = fo.V
fo.params['K_0'] = fo.K_T 
fo.params['Kdprime_0'] = -fo.params['Kprime_0']/fo.params['K_0']
fo.params['a_0'] = fo.alpha

alphas_1 = []
volumes_1 = []
for T in temperatures:
    fo.set_state(1.e5, T)
    alphas_1.append(fo.alpha)
    volumes_1.append(fo.V)
    
plt.plot(temperatures, np.array(alphas_0), label='Tref = room temperature')
plt.plot(temperatures, np.array(alphas_1), label='Tref = '+str(T_ref)+' K')

plt.title('Thermal expansion at '+str(room_pressure/1.e9)+' GPa')
plt.legend(loc='lower right')
plt.xlabel('Temperatures (K)')
plt.ylabel('Thermal expansivity (m^3/m^3/K)')
plt.show()        

plt.plot(temperatures, np.array(volumes_0), label='Tref = room temperature')
plt.plot(temperatures, np.array(volumes_1), label='Tref = '+str(T_ref)+' K')
plt.title('Thermal expansion')
plt.legend(loc='lower right')
plt.xlabel('Temperatures (K)')
plt.ylabel('Volumes (m^3/mol)')
plt.show()

pressures = np.linspace(1.e5, 50.e9, 100)
eos_1 = []
for P in pressures:
    fo.set_state(P, T_eos)
    eos_1.append(fo.V)

plt.plot(pressures/1.e9, np.array(eos_0), label='Tref = room temperature')
plt.plot(pressures/1.e9, np.array(eos_1), label='Tref = '+str(T_ref)+' K')

plt.title('Compression at '+str(T_eos)+' K')
plt.legend(loc='lower right')
plt.xlabel('Pressures (GPa)')
plt.ylabel('Volumes (m^3/mol)')
plt.show()
