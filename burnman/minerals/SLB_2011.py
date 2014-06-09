# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

SLB_2011
^^^^^^^^

Minerals from Stixrude & Lithgow-Bertelloni 2011 and references therein

"""

import burnman.mineral_helpers as bmb
from burnman.mineral import Mineral


class stishovite (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 14.02e-6,
            'K_0': 314.0e9,
            'Kprime_0': 3.8,
            'G_0': 220.0e9,
            'Gprime_0': 1.9,
            'molar_mass': .0601,
            'n': 3,
            'Debye_0': 1108.,
            'grueneisen_0': 1.37,
            'q_0': 2.8,
            'eta_s_0': 4.6}

        self.uncertainties = {
             'err_K_0':8.e9,
             'err_Kprime_0':0.1,
             'err_G_0':12.e9,
             'err_Gprime_0':0.1,
             'err_Debye_0' : 13.,
             'err_grueneisen_0': .17,
             'err_q_0': 2.2,
             'err_eta_s_0' : 1.0
            }


class periclase (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein
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
            'Debye_0': 767.,
            'grueneisen_0': 1.36,
            'q_0': 1.7, #1.7
            'eta_s_0': 2.8 } # 2.8

        self.uncertainties = {
        'err_K_0': 3.e9,
        'err_Kprime_0':.2,
        'err_G_0':1.0e9,
        'err_Gprime_0':.1,
        'err_Debye_0':9.,
        'err_grueneisen_0':.05,
        'err_q_0':.2,
        'err_eta_s_0':.2 }

class wuestite (Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 12.26e-6,
            'K_0': 179.0e9,
            'Kprime_0': 4.9,
            'G_0': 59.0e9,
            'Gprime_0': 1.4,
            'molar_mass': .0718,
            'n': 2,
            'Debye_0': 454.,
            'grueneisen_0': 1.53,
            'q_0': 1.7, #1.7
            'eta_s_0': -0.1 } #

        self.uncertainties = {
            'err_K_0':1.e9,
            'err_Kprime_0':.2,
            'err_G_0':1.e9,
            'err_Gprime_0':.1,
            'err_Debye_0':21.,
            'err_grueneisen_0':.13,
            'err_q_0':1.0,
            'err_eta_s_0':1.0}


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
    Stixrude & Lithgow-Bertelloni 2011 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 24.45e-6,
            'K_0': 251.0e9,
            'Kprime_0': 4.1,
            'G_0': 173.0e9,
            'Gprime_0': 1.7,
            'molar_mass': .1000,
            'n': 5,
            'Debye_0': 905.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'eta_s_0': 2.6 } #2.6

        self.uncertainties = {
            'err_K_0': 3.e9,
            'err_Kprime_0': 0.1,
            'err_G_0': 2.e9,
            'err_Gprime_0' : 0.0,
            'err_Debye_0': 5.,
            'err_grueneisen_0':.05,
            'err_q_0': .3,
            'err_eta_s_0':.3}


class fe_perovskite(Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 25.49e-6,
            'K_0': 272.0e9,
            'Kprime_0': 4.1,
            'G_0': 133.0e9,
            'Gprime_0': 1.4,
            'molar_mass': .1319,
            'n': 5,
            'Debye_0': 871.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'eta_s_0': 2.3 } #2.3

        self.uncertainties = {
            'err_K_0':40e9,
            'err_Kprime_0':1.,
            'err_G_0':40e9,
            'err_Gprime_0':0.0,
            'err_Debye_0':26.,
            'err_grueneisen_0':.3,
            'err_q_0':1.0,
            'err_eta_s_0':1.0}

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
