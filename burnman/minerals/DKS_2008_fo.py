# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
DKS_2008_fo
Solid and liquid forsterite from de Koker and Stixrude (2008) FPMD simulations
"""

from burnman.mineral import Mineral
from burnman.solidsolution import SolidSolution
from burnman.solutionmodel import *
from burnman.constants import gas_constant

class forsterite(Mineral):
    def __init__(self):
        K_0= 101.e9
        K_prime_0= 5.4
        K_dprime_0 = ( -143./9. - K_prime_0*(K_prime_0 - 7.)) / K_0
        n_atoms = 7.
        self.params = {
            'name': 'forsterite',
            'formula': {'Mg': 2., 'Si': 1., 'O': 4.},
            'equation_of_state': 'mg',
            'V_0': 44.8e-6, # [m^3/mol]
            'T_0': 1000., # [K]
            'E_0': 0., # [J/mol], not in article
            'S_0': 0., # [J/K/mol], nor in article
            'K_0': K_0,
            'K_prime_0': K_prime_0,
            'K_dprime_0': K_dprime_0,
            'V_x': 52.36e-6,
            'C_v_V_x': 3.1*n_atoms*gas_constant, # [J/K/mol]
            'C_v_prime_V_x': 0.05*n_atoms*gas_constant,
            'grueneisen_V_x': 1.4,
            'grueneisen_prime_V_x': 1.8
            } 
        Mineral.__init__(self)

class fo_liquid(Mineral):
    def __init__(self):
        K_0= 19.e9
        K_prime_0= 6.2
        K_dprime_0 = ( -143./9. - K_prime_0*(K_prime_0 - 7.)) / K_0
        n_atoms = 7.
        self.params = {
            'name': 'forsterite liquid',
            'formula': {'Mg': 2., 'Si': 1., 'O': 4.},
            'equation_of_state': 'mg',
            'V_0': 57.8e-6, # [m^3/mol]
            'T_0': 3000., # [K]
            'E_0': 0., # [J/mol], not in article
            'S_0': 0., # [J/K/mol], not in article
            'K_0': K_0,
            'K_prime_0': K_prime_0,
            'K_dprime_0': K_dprime_0,
            'V_x': 52.36e-6,
            'C_v_V_x': 4.29*n_atoms*gas_constant, # [J/K/mol]
            'C_v_prime_V_x': 0.68*n_atoms*gas_constant,
            'grueneisen_V_x': 0.64,
            'grueneisen_prime_V_x': -1.2
            } 
        Mineral.__init__(self)
