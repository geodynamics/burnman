# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

from minerals_base import *
                
class blank (material): #Dummy mineral so minerals.py doesn't try to import a mineral
	#list that doesn't exist. Feel free to remove this comment when minerals are added
    """
    Stixrude & Lithgow-Bertelloni 2011 and references therein 
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