# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

HP_2011
^^^^^^^^

Minerals from Holland and Powell 2011 and references therein.

Note the units in Holland and Powell's Table 1 are not SI. They are:
H_0 [kJ/mol]     i.e. multiply by 1e3 to get J/mol
err_H_0 [kJ/mol] i.e. multiply by 1e3 to get J/mol
S [J/K/mol]      i.e. S.I!
V [kJ/kbar/mol]  i.e. multiply by 1e-5 to get m^3/mol
Cp [kJ/K/mol]    is given as a + bT + cT^-2 + dT^-0.5. 
                 b is multiplied by 1e5 to be more readable. 
                 Thus, conversions to SI are [1e3, 1e-2, 1e3, 1e3]
a_0 [1e5 K^-1]   i.e. multiply by 1e-5 to get K^-1
K_0 [kbar]       i.e. multiply by 1e8 to get Pa
Kprime_0 []      i.e. S.I.!
Kdprime_0 [kbar^-1] i.e. multiply by 1e-8

The values in the database itself are NOT the same as the paper. They are strictly in the units kJ, kbar, K.

There are also parameters which deal with internal order-disorder and Landau transitions. NOTE: These still need to be implemented.
"""

from burnman.processchemistry import ProcessChemistry
from burnman.mineral import Mineral
#from burnman.solidsolution import SolidSolution

class stishovite (Mineral):
    """
    Holland and Powell, 2011 and references therein
    """
    def __init__(self):
        formula='SiO2'
        self.params = {
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -876.39e3,
            'S_0': 24.0,
            'V_0': 1.401e-5,
            'Cp': [68.1,6.01e-3,-1.9782e6,-82.1],
            'a_0': 15.8e-6,
            'K_0': 3090.e8,
            'Kprime_0': 4.6,
            'Kdprime_0': -0.00150e-8,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1]}

        self.uncertainties = {
             'err_H_0':0.49e3}

class spinel (Mineral):
    """
    Holland and Powell, 2011 and references therein
    """
    def __init__(self):
        formula='MgAl2O4'
        self.params = {
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2301.19e3,
            'S_0': 82.0,
            'V_0': 3.978e-5,
            'Cp': [222.9,6.127e-3,-1.686e6,-1551.0],
            'a_0': 1.93e-5,
            'K_0': 1922.e8,
            'Kprime_0': 4.04,
            'Kdprime_0': -0.00210e-8,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 8.00e3,
            'BW_deltaV': 0.00,
            'BW_W': 1.20e3,
            'BW_Wv': 0.00,
            'BW_n': 2.0,
            'BW_factor': 0.50}

        self.uncertainties = {
             'err_H_0':0.84e3}

class hercynite (Mineral):
    """
    Holland and Powell, 2011 and references therein
    """
    def __init__(self):
        formula='FeAl2O4'
        self.params = {
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1953.03e3,
            'S_0': 113.9,
            'V_0': 4.0750e-5,
            'Cp': [216.7,5.868e-3,-2.4302e6,-1178.3],
            'a_0': 2.06e-5,
            'K_0': 1922.e8,
            'Kprime_0': 4.04,
            'Kdprime_0': -0.00210e-8,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 18.30e3,
            'BW_deltaV': 0.00e-5,
            'BW_W': 13.60e3,
            'BW_Wv': 0.00e-5,
            'BW_n': 2.0,
            'BW_factor': 1.00}

        self.uncertainties = {
             'err_H_0':0.85e3}

class sillimanite (Mineral):
    """
    Holland and Powell, 2011 and references therein
    """
    def __init__(self):
        formula='Al2SiO5'
        self.params = {
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2585.79e3,
            'S_0': 95.4,
            'V_0': 4.986e-5,
            'Cp': [280.2,-6.9e-3,-1.3757e6,-2399.4],
            'a_0': 1.12e-5,
            'K_0': 1640.e8,
            'Kprime_0': 5.06,
            'Kdprime_0': -0.00310e-8,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 4.75e3,
            'BW_deltaV': 0.01e-5,
            'BW_W': 4.75e3,
            'BW_Wv': 0.01e-5,
            'BW_n': 1.0,
            'BW_factor': 0.25}

        self.uncertainties = {
             'err_H_0':0.68e3}

class dolomite (Mineral):
    """
    Holland and Powell, 2011 and references therein
    """
    def __init__(self):
        formula='MgCaC2O6'
        self.params = {
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -2326.22e3,
            'S_0': 156.1,
            'V_0': 6.429e-5,
            'Cp': [358.9,-4.905e-3,0.,-3456.2],
            'a_0': 3.28e-5,
            'K_0': 943.e8,
            'Kprime_0': 3.74,
            'Kdprime_0': -0.004e-8,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'BW_deltaH': 11.91e3,
            'BW_deltaV': 0.016e-5,
            'BW_W': 11.91e3,
            'BW_Wv': 0.016e-5,
            'BW_n': 1.0,
            'BW_factor': 1.0}

        self.uncertainties = {
             'err_H_0':0.58e3}


class quartz (Mineral):
    """
    Holland and Powell, 2011 and references therein
    """
    def __init__(self):
        formula='SiO2'
        self.params = {
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -910.70e3,
            'S_0': 41.43,
            'V_0': 2.269e-5,
            'Cp': [92.9,-6.42e-4,-7.149e5,-716.1],
            'a_0': 0.0e-6,
            'K_0': 730.e8,
            'Kprime_0': 6.0,
            'Kdprime_0': -0.00820e-8,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc':847.,
            'landau_Smax':4.95,
            'landau_Vmax':0.1188e-5}

        self.uncertainties = {
             'err_H_0':0.27e3}

class iron (Mineral):
    """
    Holland and Powell, 2011 and references therein
    """
    def __init__(self):
        formula='Fe'
        self.params = {
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 0.0e3,
            'S_0': 27.09,
            'V_0': 0.709e-5,
            'Cp': [46.2,5.159e-3,7.231e5,-556.2],
            'a_0': 3.56e-5,
            'K_0': 1640.e8,
            'Kprime_0': 5.16,
            'Kdprime_0': -0.00310e-8,
            'n': ProcessChemistry(formula)[0],
            'molar_mass': ProcessChemistry(formula)[1],
            'landau_Tc':1042.,
            'landau_Smax':8.3,
            'landau_Vmax':0.0000e-5}

        self.uncertainties = {
             'err_H_0':0.00e3}

