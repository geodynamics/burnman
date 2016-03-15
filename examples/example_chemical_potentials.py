# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_chemical_potentials
---------------------------

This example shows how to use the chemical potentials library of functions.

*Demonstrates:*

* How to calculate chemical potentials
* How to compute fugacities and relative fugacities

"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
import burnman.constants as constants
from burnman.processchemistry import *
from burnman.chemicalpotentials import *
import numpy as np
import burnman.minerals as minerals

if __name__ == "__main__":
    atomic_masses = read_masses()

    '''
    Here we initialise the minerals we'll be using
    '''
    P = 1.e9
    T = 1000.

    fa = minerals.HP_2011_ds62.fa()
    mt = minerals.HP_2011_ds62.mt()
    qtz = minerals.HP_2011_ds62.q()
    FMQ = [fa, mt, qtz]

    for mineral in FMQ:
        mineral.set_state(P, T)

    '''
    Here we find chemical potentials of FeO, SiO2 and O2 for
    an assemblage containing fayalite, magnetite and quartz,
    and a second assemblage of magnetite and wustite
    at 1 GPa, 1000 K
    '''

    component_formulae = ['FeO', 'SiO2', 'O2']
    component_formulae_dict = [dictionarize_formula(f)
                               for f in component_formulae]
    chem_potentials = chemical_potentials(FMQ, component_formulae_dict)

    oxygen = minerals.HP_2011_fluids.O2()
    oxygen.set_state(P, T)

    hem = minerals.HP_2011_ds62.hem()
    MH = [mt, hem]
    for mineral in MH:
        mineral.set_state(P, T)

    print('log10(fO2) at the FMQ buffer:', np.log10(fugacity(oxygen, FMQ)))
    print('log10(fO2) at the mt-hem buffer:', np.log10(fugacity(oxygen, MH)))

    print('Relative log10(fO2):', np.log10(relative_fugacity(oxygen, FMQ, MH)))

    '''
    Here we find the oxygen fugacity of the
    FMQ buffer, and compare it to published values.

    Fugacity is often defined relative to a material at
    some fixed reference pressure (in this case, O2)
    Here we use room pressure, 100 kPa
    '''

    # Set up arrays
    temperatures = np.linspace(900., 1420., 100)
    log10fO2_FMQ_ONeill1987 = np.empty_like(temperatures)
    log10fO2_FMQ = np.empty_like(temperatures)
    invT = np.empty_like(temperatures)

    # Reference and assemblage pressure
    Pr = 1.e5
    P = 1.e5
    for i, T in enumerate(temperatures):
        # Set states
        oxygen.set_state(Pr, T)
        for mineral in FMQ:
            mineral.set_state(P, T)

        # The chemical potential and fugacity of O2 at the FMQ buffer
        # according to O'Neill, 1987
        muO2_FMQ_ONeill1987 = -587474. + 1584.427 * \
            T - 203.3164 * T * np.log(T) + 0.092710 * T * T
        log10fO2_FMQ_ONeill1987[i] = np.log10(
            np.exp((muO2_FMQ_ONeill1987) / (constants.gas_constant * T)))

        invT[i] = 10000. / (T)

        # The calculated chemical potential and fugacity of O2 at the FMQ
        # buffer
        log10fO2_FMQ[i] = np.log10(fugacity(oxygen, FMQ))

    # Plot the FMQ log10(fO2) values
    plt.plot(temperatures, log10fO2_FMQ_ONeill1987,
             'k', linewidth=1., label='FMQ (O\'Neill (1987)')
    plt.plot(temperatures, log10fO2_FMQ, 'b--',
             linewidth=2., label='FMQ (HP 2011 ds62)')

    # Do the same for Re-ReO2
    '''
    Here we define two minerals, Re (rhenium) and
    ReO2 (tugarinovite)
    '''
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
    Here we find the oxygen fugacity of the Re-ReO2
    buffer, and again compare it to published values.
    '''

    # Mineral and assemblage definitions
    rhenium = Re()
    rheniumIVoxide = ReO2()
    ReReO2buffer = [rhenium, rheniumIVoxide]

    # Set up arrays
    temperatures = np.linspace(850., 1250., 100)
    log10fO2_Re_PO1994 = np.empty_like(temperatures)
    log10fO2_ReReO2buffer = np.empty_like(temperatures)

    for i, T in enumerate(temperatures):
        # Set states
        oxygen.set_state(Pr, T)
        for mineral in ReReO2buffer:
            mineral.set_state(P, T)

        # The chemical potential and fugacity of O2 at the Re-ReO2 buffer
        # according to Powncesby and O'Neill, 1994
        muO2_Re_PO1994 = -451020 + 297.595 * T - 14.6585 * T * np.log(T)
        log10fO2_Re_PO1994[i] = np.log10(
            np.exp((muO2_Re_PO1994) / (constants.gas_constant * T)))

        invT[i] = 10000. / (T)

        # The chemical potential and fugacity of O2 at the Re-ReO2 buffer
        log10fO2_ReReO2buffer[i] = np.log10(fugacity(oxygen, ReReO2buffer))

    # Plot the Re-ReO2 log10(fO2) values
    plt.plot(temperatures, log10fO2_Re_PO1994, 'k',
             linewidth=1., label='Re-ReO2 (Pownceby and O\'Neill (1994)')
    plt.plot(temperatures, log10fO2_ReReO2buffer,
             'r--', linewidth=2., label='Re-ReO2 (HP 2011 ds62)')
    plt.ylabel("log_10 (fO2)")
    plt.xlabel("T (K)")
    plt.legend(loc='lower right')
    plt.show()
