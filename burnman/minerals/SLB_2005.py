# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
SLB_2005
^^^^^^^^

Minerals from Stixrude & Lithgow-Bertelloni 2005 and references therein

"""
import burnman.mineral_helpers as bmb
from burnman.mineral import Mineral


class stishovite (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 14.02e-6,
            'K_0': 314.0e9,
            'Kprime_0': 4.4,
            'G_0': 220.0e9,
            'Gprime_0': 1.6,
            'molar_mass': .0601,
            'n': 3,
            'Debye_0': 1044.,
            'grueneisen_0': 1.34,
            'q_0': 2.4,
            'eta_s_0': 5.0 }

class periclase (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.8,
            'G_0': 131.0e9,
            'Gprime_0': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 2.8 }

class wuestite (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
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
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 24.45e-6,
            'K_0': 251.0e9,
            'Kprime_0': 4.1,
            'G_0': 175.0e9,
            'Gprime_0': 1.7,
            'molar_mass': .1000,
            'n': 5,
            'Debye_0': 1070.,
            'grueneisen_0': 1.48,
            'q_0': 1.4,
            'eta_s_0': 2.6 }

class fe_perovskite(Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 25.48e-6,
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


class mg_fe_perovskite_pt_dependent(bmb.HelperFeDependent):
    def __init__(self, iron_number_with_pt, idx):
        bmb.HelperFeDependent.__init__(self, iron_number_with_pt, idx)

    def create_inner_material(self, iron_number):
        return mg_fe_perovskite(iron_number)

class ferropericlase_pt_dependent(bmb.HelperFeDependent):
    def __init__(self, iron_number_with_pt, idx):
        bmb.HelperFeDependent.__init__(self, iron_number_with_pt, idx)

    def create_inner_material(self, iron_number):
        return ferropericlase(iron_number)
