# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


from burnman.minerals_base import *



class Speziale_fe_periclase(helper_spin_transition):  
    def __init__(self):
        helper_spin_transition.__init__(self, 60.0e9, Speziale_fe_periclase_LS(), Speziale_fe_periclase_HS())
        self.cite = 'Speziale et al. 2007'

class Speziale_fe_periclase_HS(material):
    """
    Speziale et al. 2007, Mg#=83
    """ 
    def __init__(self):
            self.params = {
                        'equation_of_state': 'mgd3',
                        'ref_V': 22.9e-6,
                        'ref_K': 157.5e9,
                        'K_prime': 3.92,
                        'molar_mass': .04567,
                        'n': 2,
                        'ref_Debye': 587,
                        'ref_grueneisen': 1.46,
                        'q0': 1.2 }

class Speziale_fe_periclase_LS(material): 
    """
    Speziale et al. 2007, Mg#=83
    """
    def __init__(self):
        self.params = {
                        'equation_of_state': 'mgd3',
                        'ref_V': 21.49e-6,
                        'ref_K': 186.0e9,
                        'K_prime': 4.6,
                        'molar_mass': .04567,
                        'n': 2,
                        'ref_Debye': 587.,
                        'ref_grueneisen': 1.46,
                        'q0': 1.2  }

    
    

        


