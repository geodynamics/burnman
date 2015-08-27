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
from scipy import optimize
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
import matplotlib.pyplot as plt

T_min_Cp = 500.
T_max_Cp = 2000.

T_max = 6000.

T_ref = 1273.
P_ref = 50.e9

T_obs = 1873.
P_obs = 50.e9


fo=minerals.HP_2011_ds62.fo()


def Cp(T, a, b, c, d):
    return a + b*T + c/T/T + d/np.sqrt(T)

fo.set_state(P_ref, T_ref)

Cps = []
temperatures = np.linspace(T_min_Cp, T_max_Cp, 100)
for T in temperatures:
    fo.set_state(P_obs, T)
    Cps.append(fo.C_p)
    
Cp_new = optimize.curve_fit(Cp, temperatures, Cps)[0]
a, b, c, d = Cp_new


fo.set_state(P_obs, T_obs)
print fo.gibbs, fo.S, fo.C_p



Cp0 = []
temperatures = np.linspace(300., T_max, 100)
for T in temperatures:
    fo.set_state(P_obs, T)
    Cp0.append(fo.C_p)


temperatures = np.linspace(300., T_max, 100)
alphas_0 = []
volumes_0 = []
gibbs_0 = []
for T in temperatures:
    fo.set_state(P_obs, T)
    alphas_0.append(fo.alpha)
    volumes_0.append(fo.V)
    gibbs_0.append(fo.gibbs)


pressures = np.linspace(1.e5, 100.e9, 100)
eos_0 = []
for P in pressures:
    fo.set_state(P, T_obs)
    eos_0.append(fo.V)


dP = 100000.
fo.set_state(P_ref-dP, T_ref)

K0 = fo.K_T 
fo.set_state(P_ref, T_ref)
K1 = fo.K_T
fo.set_state(P_ref+dP, T_ref)
K2 = fo.K_T

grad0 = (K1 - K0)/dP
grad1 = (K2 - K1)/dP


fo.set_state(P_ref, T_ref)
fo.params['T_0'] = T_ref
fo.params['P_0'] = P_ref
fo.params['H_0'] = fo.H
fo.params['S_0'] = fo.S
fo.params['V_0'] = fo.V
fo.params['Cp'] = [a, b, c, d]
fo.params['K_0'] = fo.K_T 
fo.params['a_0'] = fo.alpha
fo.params['Kprime_0'] = (K2 - K0)/(2.*dP)
fo.params['Kdprime_0'] = (grad1 - grad0)/dP


fo.set_state(P_obs, T_obs)
print fo.gibbs, fo.S, fo.C_p

Cp1 = []
temperatures = np.linspace(300., T_max, 100)
for T in temperatures:
    fo.set_state(P_obs, T)
    Cp1.append(fo.C_p)


plt.plot(temperatures, np.array(Cp0), label='pre')
plt.plot(temperatures, np.array(Cp1), label='post')
plt.legend(loc='lower right')
plt.show()


alphas_1 = []
volumes_1 = []
gibbs_1 = []
for T in temperatures:
    fo.set_state(P_obs, T)
    alphas_1.append(fo.alpha)
    volumes_1.append(fo.V)
    gibbs_1.append(fo.gibbs)
    
plt.plot(temperatures, np.array(alphas_0), label='Tref = room temperature')
plt.plot(temperatures, np.array(alphas_1), label='Tref = '+str(T_ref)+' K')

plt.title('Thermal expansion at '+str(P_obs/1.e9)+' GPa')
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

pressures = np.linspace(1.e5, 100.e9, 100)
eos_1 = []
for P in pressures:
    fo.set_state(P, T_obs)
    eos_1.append(fo.V)

plt.plot(pressures/1.e9, np.array(eos_0), label='Tref = room temperature')
plt.plot(pressures/1.e9, np.array(eos_1), label='Tref = '+str(T_ref)+' K')

plt.title('Compression at '+str(T_obs)+' K')
plt.legend(loc='lower right')
plt.xlabel('Pressures (GPa)')
plt.ylabel('Volumes (m^3/mol)')
plt.show()


plt.plot(temperatures, np.array(gibbs_0) - np.array(gibbs_1))

plt.title('Compression at '+str(T_obs)+' K')
plt.xlabel('Temperatures (K)')
plt.ylabel('diff gibbs')
plt.show()
