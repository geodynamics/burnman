# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_gibbs_modifiers
----------------

This example script demonstrates the modifications to
the gibbs free energy (and derivatives) that can be applied
as masks over the results from the equations of state.

These modifications currently take the forms:
- Landau corrections (implementations of Putnis (1992)
  and Holland and Powell (2011)
- Bragg-Williams corrections
  (implementation of Holland and Powell (1996))
- Linear (a simple delta_E + delta_V*P - delta_S*T
- Magnetic (Chin, Hertzman and Sundman (1987))

*Uses:*

* :doc:`mineral_database`


*Demonstrates:*

* creating a mineral with excess contributions
* calculating thermodynamic properties
"""
from __future__ import absolute_import

# Here we import standard python modules that are required for
# usage of BurnMan.  In particular, numpy is used for handling
# numerical arrays and mathematical operations on them, and
# matplotlib is used for generating plots of results of calculations
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))


# Here we import the relevant modules from BurnMan.  The burnman
# module imports several of the most important functionalities of
# the library, including the ability to make composites, and compute
# thermoelastic properties of them.  The minerals module includes
# the mineral physical parameters for the predefined minerals in
# BurnMan
import burnman
from burnman import minerals


if __name__ == "__main__":

    # Here we show the interesting features of Landau transitions
    # First, lets look at the P wave velocity in stishovite as it transforms
    # to the CaCl2 structure at high pressure
    stv = minerals.SLB_2011.stishovite()

    T = 1500.
    pressures = np.linspace(60.e9, 80.e9, 101)
    v_ps = np.empty_like(pressures)
    for i, P in enumerate(pressures):
        stv.set_state(P, T)
        v_ps[i] = stv.v_p

    plt.plot(pressures / 1.e9, v_ps / 1.e3, label='stishovite')
    plt.xlabel('P (GPa)')
    plt.ylabel('V_p (km/s)')
    plt.legend(loc="lower right")
    plt.show()

    # Landau transitions also cause spikes in heat capacity
    # Here we show an example of troilite, as implemented by
    # Evans et al. (2010) and incorporated into the dataset
    # of Holland and Powell (2011)

    # Here we show you how to create a mineral with a
    # Landau transition.
    # A special feature of burnman is that you can have
    # more than one Landau (or any other type of)
    # contribution.
    # Here's a copy of lot (low-temperature troilite) from
    # Holland and Powell (2011), with the Landau transition
    # of tro also included.
    from burnman.processchemistry import read_masses, dictionarize_formula, formula_mass
    atomic_masses = read_masses()

    class lot (burnman.Mineral):

        def __init__(self):
            formula = 'Fe1.0S1.0'
            formula = dictionarize_formula(formula)
            self.params = {
                'name': 'lot',
                'formula': formula,
                'equation_of_state': 'hp_tmt',
                'H_0': -102160.0,
                'S_0': 60.0,
                'V_0': 1.818e-05,
                'Cp': [50.2, 0.011052, -940000.0, 0.0],
                'a_0': 4.93e-05,
                'K_0': 65800000000.0,
                'Kprime_0': 4.17,
                'Kdprime_0': -6.3e-11,
                'n': sum(formula.values()),
                'molar_mass': formula_mass(formula, atomic_masses)}
            self.property_modifiers = [
                ['landau_hp', {'P_0': 100000.0,
                               'T_0': 298.15,
                               'Tc_0': 420.0,
                               'S_D': 10.0,
                               'V_D': 0.0}],
                ['landau_hp', {'P_0': 100000.0,
                               'T_0': 298.15,
                               'Tc_0': 598.0,
                               'S_D': 12.0,
                               'V_D': 4.1e-7}]]
            burnman.Mineral.__init__(self)

    troilite = lot()
    lot = minerals.HP_2011_ds62.lot()
    tro = minerals.HP_2011_ds62.tro()
    P = 1.e5
    temperatures = np.linspace(300., 1300., 101)
    C_ps_troilite = np.empty_like(temperatures)
    C_ps_lot = np.empty_like(temperatures)
    C_ps_tro = np.empty_like(temperatures)
    for i, T in enumerate(temperatures):
        troilite.set_state(P, T)
        C_ps_troilite[i] = troilite.C_p
        lot.set_state(P, T)
        C_ps_lot[i] = lot.C_p
        tro.set_state(P, T)
        C_ps_tro[i] = tro.C_p

    plt.plot(temperatures, C_ps_lot, 'r--', label='low temperature (HP2011)')
    plt.plot(temperatures, C_ps_tro, 'g--', label='high temperature (HP2011)')
    plt.plot(temperatures, C_ps_troilite, 'b-', label='troilite')
    plt.xlabel('T (K)')
    plt.ylabel('C_p (J/K/mol)')
    plt.legend(loc="lower right")
    plt.show()

    # Spinel is a mineral with a Bragg-Williams type model
    sp = minerals.HP_2011_ds62.sp()
    P = 1.e5
    temperatures = np.linspace(300., 1300., 101)
    C_ps = np.empty_like(temperatures)
    for i, T in enumerate(temperatures):
        sp.set_state(P, T)
        C_ps[i] = sp.C_p
        # print sp._property_modifiers

    plt.plot(temperatures, C_ps, label='spinel')
    plt.xlabel('T (K)')
    plt.ylabel('C_p (J/K/mol)')
    plt.legend(loc="lower right")
    plt.show()

    # Wuestite has a Landau-type transition at low temperature,
    # but we could also choose to simplify things by just having an excess entropy
    # to estimate the thermal properties at high temperature
    # Here we ignore the 0 Pa, 0 K gibbs and volume contributions, as the endmember
    # properties would need refitting too...
    class wuestite (burnman.Mineral):

        def __init__(self):
            formula = 'FeO'
            formula = dictionarize_formula(formula)
            self.params = {
                'name': 'Wuestite',
                'formula': formula,
                'equation_of_state': 'slb3',
                'F_0': -242000.0,
                'V_0': 1.226e-05,
                'K_0': 1.79e+11,
                'Kprime_0': 4.9,
                'Debye_0': 454.0,
                'grueneisen_0': 1.53,
                'q_0': 1.7,
                'G_0': 59000000000.0,
                'Gprime_0': 1.4,
                'eta_s_0': -0.1,
                'n': sum(formula.values()),
                'molar_mass': formula_mass(formula, atomic_masses)}

            self.property_modifiers = [
                ['linear', {'delta_E': 0., 'delta_S': 12., 'delta_V': 0.}]]

            self.uncertainties = {
                'err_F_0': 1000.0,
                'err_V_0': 0.0,
                'err_K_0': 1000000000.0,
                'err_K_prime_0': 0.2,
                'err_Debye_0': 21.0,
                'err_grueneisen_0': 0.13,
                'err_q_0': 1.0,
                'err_G_0': 1000000000.0,
                'err_Gprime_0': 0.1,
                'err_eta_s_0': 1.0}
            burnman.Mineral.__init__(self)

    wus = wuestite()
    wus_HP = burnman.minerals.HP_2011_ds62.fper()

    P = 1.e5
    temperatures = np.linspace(300., 1300., 101)
    Ss = np.empty_like(temperatures)
    Ss_HP = np.empty_like(temperatures)
    for i, T in enumerate(temperatures):
        wus.set_state(P, T)
        Ss[i] = wus.S
        wus_HP.set_state(P, T)
        Ss_HP[i] = wus_HP.S

    plt.plot(temperatures, Ss, label='linear')
    plt.plot(temperatures, Ss_HP, label='HP_2011_ds62')
    plt.xlabel('T (K)')
    plt.ylabel('S (J/K/mol)')
    plt.legend(loc="lower right")
    plt.show()
