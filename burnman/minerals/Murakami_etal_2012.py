# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

from burnman.minerals_base import *






class mg_perovskite(material):  
    """
    Murakami et al. (2012) supplementary table 5 and references therein, V_0 from Stixrude & Lithgow-Bertolloni 2005
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 24.45e-6,  # S & L-B 2005
            'K_0': 281e9,
            'Kprime_0': 4.1,
            'G_0': 173e9,
            'Gprime_0': 1.56,
            'molar_mass': .100,
            'n': 5,
            'Debye_0': 1070.,
            'grueneisen_0': 1.48,
            'q_0': 1.4,
            'eta_s_0': 2.4 } 

class mg_perovskite_3rdorder(material):
    """
    Murakami et al. (2012) third order fit to supplementary table 1, includes 4% Al 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 24.45e-6,  # S & L-B 2005
            'K_0': 281e9,
            'Kprime_0': 4.1,
            'G_0': 171.42e9,
            'Gprime_0': 1.83,
            'molar_mass': .100,
            'n': 5,
            'Debye_0': 1070.,
            'grueneisen_0': 1.48,
            'q_0': 1.4,
            'eta_s_0': 2.4 }

class fe_perovskite(material): 
    """
    Murakami et al. (2012), personal communication, Mg#=94, Al=4%
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 24.607e-6,
            'K_0': 251.9e9, 
            'Kprime_0': 4.01,
            'G_0': 164.7e9,
            'Gprime_0': 1.58,
            'molar_mass': .102,
            'n': 5,
            'Debye_0': 1054.,
            'grueneisen_0': 1.48,
            'q_0': 1.4,
            'eta_s_0': 2.4 } 


            
class mg_periclase(material):
    """
    Murakami et al. (2012) supplementary table 5 and references therein, V_0 from Stixrude & Lithgow-Bertolloni 2005
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 11.24e-6,  # S & L-B 2005
            'K_0': 161e9,
            'Kprime_0': 3.9,
            'G_0': 131e9,
            'Gprime_0': 1.92,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773., # S& L-B 2005
            'grueneisen_0': 1.5,
            'q_0': 1.5, #S&L-B 2005    
            'eta_s_0': 3.0 }             



class fe_periclase(helper_spin_transition):
    def __init__(self):
        helper_spin_transition.__init__(self, 63.0e9, fe_periclase_LS(), fe_periclase_HS())

class fe_periclase_3rd(helper_spin_transition):
    def __init__(self):
        helper_spin_transition.__init__(self, 63.0e9, fe_periclase_LS(), fe_periclase_HS())

class fe_periclase_HS(material):  # From Murakami's emails, see Cayman for details, represents Mg# = .79
    """
    Murakami et al. (2012), personal communication, Mg#=79 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 11.412e-6,
            'K_0': 159.1e9,
            'Kprime_0': 4.11,
            'G_0': 105.43e9,
            'Gprime_0': 1.773,
            'molar_mass': .047,
            'n': 2,
            'Debye_0': 706.,
            'grueneisen_0': 1.45,
            'q_0': 1.5,
            'eta_s_0': 2.54 }

class fe_periclase_LS(material):  # From Murakami's emails, see Cayman for details, represents Mg# = .79
    """
    Murakami et al. (2012), personal communication, Mg#=79
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 11.171e-6,
            'K_0': 170.0e9,
            'Kprime_0': 4.00,
            'G_0': 116.34e9,
            'Gprime_0': 1.668,
            'molar_mass': .047,
            'n': 2,
            'Debye_0': 706.,
            'grueneisen_0': 1.45,
            'q_0': 1.5,
            'eta_s_0': 2.54}



class fe_periclase_HS_3rd(material): 
    """
    Murakami et al. (2012), personal communication, Mg#=92
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 11.412e-6,
            'K_0': 159.1e9,
            'Kprime_0': 4.11,
            'G_0': 129.35e9,
            'Gprime_0': 1.993,
            'molar_mass': .0469,
            'n': 2,
            'Debye_0': 706.,
            'grueneisen_0': 1.45,
            'q_0': 1.5,
            'eta_s_0': 2.54 }

class fe_periclase_LS_3rd(material):  
    """
    Murakami et al. (2012), personal communication, Mg#=92
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 11.171e-6,
            'K_0': 170.0e9,
            'Kprime_0': 4,
            'G_0': 151.67e9,
            'Gprime_0': 1.754,
            'molar_mass': .0469,
            'n': 2,
            'Debye_0': 706.,
            'grueneisen_0': 1.45,
            'q_0': 1.5,
            'eta_s_0': 2.54}
                
