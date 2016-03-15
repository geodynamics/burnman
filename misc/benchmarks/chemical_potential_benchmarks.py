from __future__ import absolute_import
# Benchmarks for the chemical potential functions
import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))

import numpy as np
import matplotlib.pyplot as plt

import burnman
from burnman.processchemistry import *
from burnman.chemicalpotentials import *
import burnman.constants as constants
atomic_masses = read_masses()


class Re (burnman.Mineral):

    def __init__(self):
        formula = 'Re1.0'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Re',
            'formula': formula,
            'equation_of_state': 'hp_tmt',
            'H_0': 0.0,
            'S_0': 36.53,
            'V_0': 8.862e-06,
            'Cp': [23.7, 0.005448, 68.0, 0.0],
            'a_0': 1.9e-05,
            'K_0': 3.6e+11,
            'Kprime_0': 4.05,
            'Kdprime_0': -1.1e-11,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}
        burnman.Mineral.__init__(self)


class ReO2 (burnman.Mineral):

    def __init__(self):
        formula = 'Re1.0O2.0'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'ReO2',
            'formula': formula,
            'equation_of_state': 'hp_tmt',
            'H_0': -445140.0,
            'S_0': 47.82,
            'V_0': 1.8779e-05,
            'Cp': [76.89, 0.00993, -1207130.0, -208.0],
            'a_0': 4.4e-05,
            'K_0': 1.8e+11,
            'Kprime_0': 4.05,
            'Kdprime_0': -2.25e-11,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}
        burnman.Mineral.__init__(self)

'''
Here we find the oxygen fugacity of the FMQ assemblage
and also the Re-ReO2 buffer

Fugacity is often defined relative to a material at
some fixed reference pressure
Here we use room pressure, 100 kPa
'''
fa = burnman.minerals.HP_2011_ds62.fa()
mt = burnman.minerals.HP_2011_ds62.mt()
qtz = burnman.minerals.HP_2011_ds62.q()
FMQ = [fa, mt, qtz]

oxygen = burnman.minerals.HP_2011_fluids.O2()

rhenium = Re()
rheniumIVoxide = ReO2()
ReReO2buffer = [rhenium, rheniumIVoxide]

Pr = 1.e5

temperatures = np.linspace(900., 1420., 100)
log10fO2_FMQ_ONeill1987 = np.empty_like(temperatures)
log10fO2_FMQ = np.empty_like(temperatures)
invT = np.empty_like(temperatures)

P = 1.e5
for i, T in enumerate(temperatures):
    oxygen.set_state(Pr, T)
    for mineral in FMQ:
        mineral.set_state(P, T)
    for mineral in ReReO2buffer:
        mineral.set_state(P, T)

    muO2_FMQ_ONeill1987 = -587474. + 1584.427 * \
        T - 203.3164 * T * np.log(T) + 0.092710 * T * T
    log10fO2_FMQ_ONeill1987[i] = np.log10(
        np.exp((muO2_FMQ_ONeill1987) / (constants.gas_constant * T)))

    invT[i] = 10000. / (T)
    log10fO2_FMQ[i] = np.log10(fugacity(oxygen, FMQ))

plt.plot(temperatures, log10fO2_FMQ_ONeill1987,
         'k', linewidth=3., label='FMQ (O\'Neill (1987)')
plt.plot(temperatures, log10fO2_FMQ, 'b--',
         linewidth=3., label='FMQ (HP 2011 ds62)')


temperatures = np.linspace(850., 1250., 100)
log10fO2_Re_PO1994 = np.empty_like(temperatures)
log10fO2_ReReO2buffer = np.empty_like(temperatures)

for i, T in enumerate(temperatures):

    oxygen.set_state(Pr, T)
    for mineral in FMQ:
        mineral.set_state(P, T)
    for mineral in ReReO2buffer:
        mineral.set_state(P, T)

    muO2_Re_PO1994 = -451020 + 297.595 * T - 14.6585 * T * np.log(T)
    log10fO2_Re_PO1994[i] = np.log10(
        np.exp((muO2_Re_PO1994) / (constants.gas_constant * T)))

    invT[i] = 10000. / (T)
    log10fO2_ReReO2buffer[i] = np.log10(fugacity(oxygen, ReReO2buffer))

plt.plot(temperatures, log10fO2_Re_PO1994, 'k',
         linewidth=3., label='Re-ReO2 (Pownceby and O\'Neill (1994)')
plt.plot(temperatures, log10fO2_ReReO2buffer, 'r--',
         linewidth=3., label='Re-ReO2 (HP 2011 ds62)')
plt.ylabel("log_10 (fO2)")
plt.xlabel("T (K)")
plt.legend(loc='lower right')
plt.show()
