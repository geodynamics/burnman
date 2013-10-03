# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

from burnman.minerals_base import *
                




class stishovite (helper_uncertainty):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self,peturb=[0.0,0.0,0.0,0.0]):
        self.params = {
            'equation_of_state': 'slb3',
            'ref_V': 14.02e-6,
            'ref_K': 314.0e9,
            'err_ref_K':8.e9,
            'K_prime': 3.8,
            'err_K_prime':0.1,
            'ref_G': 220.0e9,
            'err_ref_G':12.e9,
            'G_prime': 1.9,
            'err_G_prime':0.1,
            'molar_mass': .0601,
            'n': 3,
            'ref_Debye': 1108.,
            'err_ref_Debye' : 13.,
            'ref_grueneisen': 1.37,
            'err_ref_grueneisen': .17,
            'q0': 2.8,
            'err_q0': 1.2,    # decreased so things don't crash... not the published value (which is 2.2)
            'eta_0s': 4.6,
            'err_eta_0s' : 1.0 }
        helper_uncertainty.__init__(self,peturb)


class periclase (material):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 11.24e-6,
            'ref_K': 161.0e9,
            'K_prime': 3.8,
            'ref_G': 131.0e9,
            'G_prime': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'ref_Debye': 767.,
            'ref_grueneisen': 1.36,
            'q0': 1.7,
            'eta_0s': 2.8 }

class wuestite (material):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 12.26e-6,
            'ref_K': 179.0e9,
            'K_prime': 4.9,
            'ref_G': 59.0e9,
            'G_prime': 1.4,
            'molar_mass': .0718,
            'n': 2,
            'ref_Debye': 454.,
            'ref_grueneisen': 1.53,
            'q0': 1.7,
            'eta_0s': 1.4 }



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
            'ref_V': 24.45e-6,
            'ref_K': 251.0e9,   
            'K_prime': 4.1,     
            'ref_G': 173.0e9,  
            'G_prime': 1.7,  
            'molar_mass': .1000,
            'n': 5,
            'ref_Debye': 905.,
            'ref_grueneisen': 1.57,
            'q0': 1.1,
            'eta_0s': 2.6 }

class fe_perovskite(material):
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 25.49e-6,
            'ref_K': 272.0e9, 
            'K_prime': 4.1,  
            'ref_G': 133.0e9,
            'G_prime': 1.4,   
            'molar_mass': .1319, 
            'n': 5,
            'ref_Debye': 871.,
            'ref_grueneisen': 1.57,
            'q0': 1.1,
            'eta_0s': 2.3 }


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
