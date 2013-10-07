# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

from burnman.minerals_base import *
                




class stishovite (helper_uncertainty):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self,peturb=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 14.02e-6,
            'K_0': 314.0e9,
            'err_K_0':8.e9,
            'Kprime_0': 3.8,
            'err_Kprime_0':0.1,
            'G_0': 220.0e9,
            'err_G_0':12.e9,
            'Gprime_0': 1.9,
            'err_Gprime_0':0.1,
            'molar_mass': .0601,
            'n': 3,
            'Debye_0': 1108.,
            'err_Debye_0' : 13.,
            'grueneisen_0': 1.37,
            'err_grueneisen_0': .17,
            'q_0': 2.8,
            'err_q_0': 1.2,    # decreased so things don't crash... not the published value (which is 2.2)
            'eta_s_0': 4.6,
            'err_eta_s_0' : 1.0 }
        helper_uncertainty.__init__(self,peturb)


class periclase (material):
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
            'q_0': 1.7,
            'eta_s_0': 2.8 }

class wuestite (material):
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
            'q_0': 1.7,
            'eta_s_0': 1.4 }



class ferropericlase(helper_solid_solution):
    def __init__(self, fe_num):
        base_materials = [periclase(), wuestite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        helper_solid_solution.__init__(self, base_materials, molar_fraction)


class mg_fe_perovskite(helper_solid_solution):
    def __init__(self, fe_num):
        base_materials = [mg_perovskite(), fe_perovskite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        helper_solid_solution.__init__(self, base_materials, molar_fraction)

class mg_perovskite(material):
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
            'eta_s_0': 2.6 }

class fe_perovskite(material):
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
            'eta_s_0': 2.3 }


class mg_fe_perovskite_pt_dependent(helper_fe_dependent):
    def __init__(self, iron_number_with_pt, idx):
        helper_fe_dependent.__init__(self, iron_number_with_pt, idx)

    def create_inner_material(self, iron_number):
        return mg_fe_perovskite(iron_number)

class ferropericlase_pt_dependent(helper_fe_dependent):
    def __init__(self, iron_number_with_pt, idx):
        helper_fe_dependent.__init__(self, iron_number_with_pt, idx)

    def create_inner_material(self, iron_number):
        return ferropericlase(iron_number)
