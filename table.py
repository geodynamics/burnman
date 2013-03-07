# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

""" Generates a text table with mineral properties. """

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals
from burnman import tools

if __name__ == "__main__":    
    
    phases = [minerals.stishovite(), \
    minerals.periclase(), \
    minerals.wustite(), \
    minerals.ferropericlase(0), \
    minerals.mg_fe_perovskite(0), \
    minerals.mg_perovskite(), \
    minerals.fe_perovskite(), \
    minerals.Speziale_fe_periclase_LS(), \
    minerals.Speziale_fe_periclase_HS(), \
    minerals.Catalli_fe_perovskite_LS(), \
    minerals.Catalli_fe_perovskite_HS(), \
    minerals.Murakami_fe_perovskite(), \
    minerals.Murakami_fe_periclase_LS(), \
    minerals.Murakami_fe_periclase_HS()   ]
    #mg_fe_perovskite_pt_dependent
    #ferropericlase_pt_dependent
    
    params = ['ref_V','ref_K','K_prime','ref_mu','mu_prime','molar_mass','n','ref_Debye','ref_grueneisen','q0','eta_0s']
    
    
    def create_list(mineral):
        row = [ p.__class__.__name__ ]
        for param in params:
            
            row.append(str(p.params[param] if param in p.params else ""))
        return row
    
    
    
    table = [ ['Name'] + params ]
    
    for p in phases:
        p.set_method('bm')
        p.set_state(0,300)
        row = create_list(p)
        table.append(row)
            
    
    tools.pretty_print_table(table, False)
    
