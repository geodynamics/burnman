# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


from minerals_base import *



class mg_perovskite(material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 24.43e-6,
            'ref_K': 250.0e9,   
            'K_prime': 4.0,     
            'ref_G': 175.0e9,  
            'G_prime': 1.8,    
            'molar_mass': .1020,
            'n': 5,
            'ref_Debye': 1070.,
            'ref_grueneisen': 1.48,
            'q0': 1.4} 

class fe_perovskite(material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 25.34e-6,
            'ref_K': 250.0e9, 
            'K_prime': 4.0,  
            'ref_G': 135.0e9, 
            'G_prime': 1.3,  
            'molar_mass': .1319, 
            'n': 5,
            'ref_Debye': 841.,
            'ref_grueneisen': 1.48,
            'q0': 1.4} 

class al_perovskite(material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 24.58e-6,
            'ref_K': 249.0e9,
            'K_prime': 4.0,
            'ref_G': 165.0e9,
            'G_prime': 1.8,
            'molar_mass': .1005,
            'n': 5,
            'ref_Debye': 1021.,
            'ref_grueneisen': 1.48,
            'q0': 1.4}

class ca_perovskite(material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 27.45e-6,
            'ref_K': 236.0e9,
            'K_prime': 3.9,
            'ref_G': 165.0e9,
            'G_prime': 2.46,
            'molar_mass': .11616,
            'n': 5,
            'ref_Debye': 984.,
            'ref_grueneisen': 1.53,
            'q0': 1.6}


class periclase (material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 11.25e-6,
            'ref_K': 160.1e9,
            'K_prime': 3.83,
            'ref_G': 130.0e9,
            'G_prime': 2.2,
            'molar_mass': .0403,
            'n': 2,
            'ref_Debye': 673.,
            'ref_grueneisen': 1.41,
            'q0': 1.3 }

class wuestite (material): # Matas et al 2007 Tables 1&2
    """
    Matas et al. 2007 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'mgd2',
            'ref_V': 12.26e-6,
            'ref_K': 160.1e9,
            'K_prime': 3.83,
            'ref_G': 46.0e9,
            'G_prime':  0.6,
            'molar_mass': .0718,
            'n': 2,
            'ref_Debye': 673.,
            'ref_grueneisen': 1.41,
            'q0': 1.3 }


