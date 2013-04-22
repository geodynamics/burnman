# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

from minerals_base import *






class Murakami_mg_perovskite(material):  
    """
    Murakami et al. (2012) supplementary table 5 and references therein, ref_V from Stixrude & Lithgow-Bertolloni 2005
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 24.45e-6,  # S & L-B 2005
            'ref_K': 281e9,
            'K_prime': 4.1,
            'ref_mu': 173e9,
            'mu_prime': 1.56,
            'molar_mass': .100,
            'n': 5,
            'ref_Debye': 1070.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.4 } 

class Murakami_mg_perovskite_3rdorder(material):
    """
    Murakami et al. (2012) third order fit to supplementary table 1, includes 4% Al 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 24.45e-6,  # S & L-B 2005
            'ref_K': 281e9,
            'K_prime': 4.1,
            'ref_mu': 171.42e9,
            'mu_prime': 1.83,
            'molar_mass': .100,
            'n': 5,
            'ref_Debye': 1070.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 2.4 }

class Murakami_fe_perovskite(material): 
    """
    Murakami et al. (2012), personal communication, Mg#=94, Al=4%
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 24.607e-6,
            'ref_K': 251.9e9, 
            'K_prime': 4.01,
            'ref_mu': 164.7e9,
            'mu_prime': 1.58,
            'molar_mass': .1053,
            'n': 5,
            'ref_Debye': 1054.,
            'ref_grueneisen': 1.48,
            'q0': 1.4,
            'eta_0s': 1.48 } 


            
class Murakami_mg_periclase(material):
    """
    Murakami et al. (2012) supplementary table 5 and references therein, ref_V from Stixrude & Lithgow-Bertolloni 2005
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 11.24e-6,  # S & L-B 2005
            'ref_K': 161e9,
            'K_prime': 3.9,
            'ref_mu': 131e9,
            'mu_prime': 1.92,
            'molar_mass': .0403,
            'n': 2,
            'ref_Debye': 773., # S& L-B 2005
            'ref_grueneisen': 1.5,
            'q0': 1.5, #S&L-B 2005    
            'eta_0s': 3.0 }             



class Murakami_fe_periclase(helper_spin_transition):
    def __init__(self):
        helper_spin_transition.__init__(self, 63.0e9, Murakami_fe_periclase_LS(), Murakami_fe_periclase_HS())

class Murakami_fe_periclase_3rd(helper_spin_transition):
    def __init__(self):
        helper_spin_transition.__init__(self, 63.0e9, Murakami_fe_periclase_LS(), Murakami_fe_periclase_HS())

class Murakami_fe_periclase_HS(material):  # From Murakami's emails, see Cayman for details, represents Mg# = .79
    """
    Murakami et al. (2012), personal communication, Mg#=79 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 11.412e-6,
            'ref_K': 159.1e9,
            'K_prime': 4.11,
            'ref_mu': 105.43e9,
            'mu_prime': 1.773,
            'molar_mass': .0469,
            'n': 2,
            'ref_Debye': 706.,
            'ref_grueneisen': 1.45,
            'q0': 1.5,
            'eta_0s': 2.54 }

class Murakami_fe_periclase_LS(material):  # From Murakami's emails, see Cayman for details, represents Mg# = .79
    """
    Murakami et al. (2012), personal communication, Mg#=79
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'ref_V': 11.171e-6,
            'ref_K': 170.0e9,
            'K_prime': 4,
            'ref_mu': 116.34e9,
            'mu_prime': 1.668,
            'molar_mass': .0469,
            'n': 2,
            'ref_Debye': 706.,
            'ref_grueneisen': 1.45,
            'q0': 1.5,
            'eta_0s': 2.54}



class Murakami_fe_periclase_HS_3rd(material): 
    """
    Murakami et al. (2012), personal communication, Mg#=92
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 11.412e-6,
            'ref_K': 159.1e9,
            'K_prime': 4.11,
            'ref_mu': 129.35e9,
            'mu_prime': 1.993,
            'molar_mass': .0469,
            'n': 2,
            'ref_Debye': 706.,
            'ref_grueneisen': 1.45,
            'q0': 1.5,
            'eta_0s': 2.54 }

class Murakami_fe_periclase_LS_3rd(material):  
    """
    Murakami et al. (2012), personal communication, Mg#=92
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 11.171e-6,
            'ref_K': 170.0e9,
            'K_prime': 4,
            'ref_mu': 151.67e9,
            'mu_prime': 1.754,
            'molar_mass': .0469,
            'n': 2,
            'ref_Debye': 706.,
            'ref_grueneisen': 1.45,
            'q0': 1.5,
            'eta_0s': 2.54}
                
