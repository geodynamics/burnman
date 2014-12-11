import os, sys, numpy as np, matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
from chemicalpotentials import *
import numpy as np

'''
Here we initialise the minerals we'll be using
'''
P=1.e8
T=1000.

fa=minerals.HP_2011_ds62.fa()
mt=minerals.HP_2011_ds62.mt()
qtz=minerals.HP_2011_ds62.q()
FMQ=[fa, mt, qtz]

for mineral in FMQ:
    mineral.set_method('mtait')
    mineral.set_state(P, T)


'''
Here we find chemical potentials of FeO, SiO2 and O2 for
an assemblage containing fayalite, magnetite and quartz
at 0.1 GPa, 1000 K
'''


component_formulae=['FeO', 'SiO2', 'O2']
component_formulae_dict=[dictionarize_formula(f) for f in component_formulae]
chem_potentials=chemicalpotentials(FMQ, component_formulae_dict)

for idx, component_formula in enumerate(component_formulae):
    print "mu("+component_formula+"):", chem_potentials[idx]
print ''


'''
Here we find the oxygen fugacity of the FMQ assemblage
and also the Re-ReO2 buffer

Fugacity is often defined relative to a material at 
some fixed reference pressure
Here we use room pressure, 100 kPa
'''
O2_gas=minerals.HP_2011_fluids.O2()
O2_gas.set_method('cork')

rhenium=minerals.Metal_Metal_oxides.Re()
rheniumIVoxide=minerals.Metal_Metal_oxides.ReO2()
ReReO2buffer=[rhenium, rheniumIVoxide]

for mineral in ReReO2buffer:
    mineral.set_method('mtait')


Pr=1.e5
O2=dictionarize_formula('O2')

temperatures = np.linspace(200, 2400, 100)
log10fO2_FMQ = np.empty_like(temperatures)
log10fO2_ReReO2buffer = np.empty_like(temperatures)
invT = np.empty_like(temperatures)
for i, T in enumerate(temperatures):
    O2_gas.set_state(Pr, T+273.15)
    for mineral in FMQ:
        mineral.set_state(P, T+273.15)
    for mineral in ReReO2buffer:
        mineral.set_state(P, T+273.15)
    invT[i] = 10000./(T+273.15)
    log10fO2_FMQ[i] = np.log10(fugacity(O2, O2_gas, FMQ))
    log10fO2_ReReO2buffer[i] = np.log10(fugacity(O2, O2_gas, ReReO2buffer))

plt.plot(invT, log10fO2_FMQ, 'b--', linewidth=3.)
plt.plot(invT, log10fO2_ReReO2buffer, 'r--', linewidth=3.)
plt.xlim(4.0,20.0)
plt.ylim(-36.,4.)
plt.ylabel("log_10 (fO2)")
plt.xlabel("1e4/T (K^-1)")
plt.show()
