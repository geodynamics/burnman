# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

from minerals_base import *
                





class stishovite (material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'ref_V': 14.02e-6,
            'ref_K': 314.0e9,
            'K_prime': 4.4,
            'ref_mu': 220.0e9,
            'mu_prime': 1.6,
            'molar_mass': .0601,
            'n': 3,
            'ref_Debye': 1044.,
            'ref_grueneisen': 1.34,
            'q0': 2.4,
            'eta_0s': 5.0 }

class periclase (material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 11.24e-6,
            'ref_K': 161.0e9,
            'K_prime': 3.8,
            'ref_mu': 131.0e9,
            'mu_prime': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'ref_Debye': 773.,
            'ref_grueneisen': 1.5,
            'q0': 1.5,
            'eta_0s': 2.8 }

class wuestite (material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 12.06e-6,
            'ref_K': 152.0e9,
            'K_prime': 4.9,
            'ref_mu': 47.0e9,
            'mu_prime': 0.7,
            'molar_mass': .0718,
            'n': 2,
            'ref_Debye': 455.,
            'ref_grueneisen': 1.28,
            'q0': 1.5,
            'eta_0s': 0.8 }



class ferropericlase(helper_volumetric_mixing):
    def __init__(self, fe_num):
        base_materials = [periclase(), wuestite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        helper_volumetric_mixing.__init__(self, base_materials, molar_fraction)



class mg_fe_perovskite(helper_volumetric_mixing):
    def __init__(self, fe_num):
        base_materials = [mg_perovskite(), fe_perovskite()]
        molar_fraction = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        helper_volumetric_mixing.__init__(self, base_materials, molar_fraction)


class mg_perovskite(material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein  
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 24.45e-6,
            'ref_K': 251.0e9,   
            'K_prime': 4.1,     
            'ref_mu': 175.0e9,  
            'mu_prime': 1.7,  
            'molar_mass': .1000,
            'n': 5,
            'ref_Debye': 1070.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.6 }

class fe_perovskite(material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 25.48e-6,
            'ref_K': 281.0e9, 
            'K_prime': 4.1,  
            'ref_mu': 138.0e9,
            'mu_prime': 1.7,   
            'molar_mass': .1319, 
            'n': 5,
            'ref_Debye': 841.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.1 }


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
