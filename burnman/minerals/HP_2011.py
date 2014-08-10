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
S [J/mol]        i.e. S.I!
V [kJ/kbar/mol]  i.e. multiply by 1e-5 to get m^3/mol
Cp [kJ/K/mol]    is given as a + bT + cT^-2 + dT^-0.5. 
                 b is multiplied by 1e5 to be more readable. 
                 Thus, conversions to SI are [1e3, 1e-2, 1e3, 1e3]
a_0 [1e5 K^-1]   i.e. multiply by 1e-5 to get K^-1
K_0 [kbar]       i.e. multiply by 1e8 to get Pa
Kprime_0 []      i.e. SI!
Kdprime_0 [kbar^-1] i.e. multiply by 1e-8

The values in the database itself are NOT the same as the paper. They are strictly in the units kJ, kbar, K.

There are also parameters which deal with internal order-disorder and Landau transitions. NOTE: These still need to be implemented.
"""

from burnman.mineral import Mineral
from burnman.solidsolution import SolidSolution

class stishovite (Mineral):
    """
    Holland and Powell, 2011 and references therein
    """
    def __init__(self):
        self.params = {
            'formula': 'SiO2',
            'equation_of_state': 'mtait',
            'H_0': -876.39e3,
            'S_0': 24.0,
            'V_0': 1.401e-5,
            'Cp': [68.1,6.01e-3,-1.9782e6,-82.1],
            'a_0': 15.8e-6
            'K_0': 3090.e8
            'Kprime_0': 4.6,
            'Kdprime_0': -0.00150e-8,
            'molar_mass': .0601,
            'n': 3}

        self.uncertainties = {
             'err_H_0':0.49e3}

