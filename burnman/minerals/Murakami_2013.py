# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
Murakami_2013
^^^^^^^^^^^^^

Minerals from Murakami 2013 and references therein

"""

import burnman.mineral_helpers as bmb
from burnman.mineral import Mineral





class periclase (Mineral):
    """
    Murakami 2013 and references therine
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.9,
            'G_0': 130.9e9,
            'Gprime_0': 1.92,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 2.3 }

class wuestite (Mineral):
    """
    Muarakami 2013 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 12.06e-6,
            'K_0': 152.0e9,
            'Kprime_0': 4.9,
            'G_0': 47.0e9,
            'Gprime_0': 0.7,
            'molar_mass': .0718,
            'n': 2,
            'Debye_0': 455.,
            'grueneisen_0': 1.28,
            'q_0': 1.5,
            'eta_s_0': 0.8 }



class ferropericlase(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        base_materials = [periclase(), wuestite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)



class mg_fe_perovskite(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        base_materials = [mg_perovskite(), fe_perovskite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)


class mg_perovskite(Mineral):
    """
    Murakami 2013 and references therin
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 24.45e-6,
            'K_0': 253.0e9,
            'Kprime_0': 4.1,
            'G_0': 172.9e9,
            'Gprime_0': 1.56,
            'molar_mass': .1000,
            'n': 5,
            'Debye_0': 1100.,
            'grueneisen_0': 1.4,
            'q_0': 1.4,
            'eta_s_0': 2.6 }

class fe_perovskite(Mineral):
    """
    Murakami 2013 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 25.49e-6,
            'K_0': 281.0e9,
            'Kprime_0': 4.1,
            'G_0': 138.0e9,
            'Gprime_0': 1.7,
            'molar_mass': .1319,
            'n': 5,
            'Debye_0': 841.,
            'grueneisen_0': 1.48,
            'q_0': 1.4,
            'eta_s_0': 2.1 }


