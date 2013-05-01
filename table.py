# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

""" Generates a text table with mineral properties. """

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
import inspect

from burnman import minerals
from burnman import tools

if __name__ == "__main__":    

    names = set()
    phases = []

    libs = dir(minerals)
    for l in libs:
        mineralgroup = getattr(minerals, l)
        
        if mineralgroup.__class__.__name__ == "module":
            for m in dir(mineralgroup):
                mineral = getattr(mineralgroup, m)
                #print mineral
                if inspect.isclass(mineral) and mineral!=burnman.material and issubclass(mineral, burnman.material):
                    #print mineral.__module__ + mineral.__name__
                    name = mineral.__module__.replace("minlib_","").replace("burnman.","minerals.") + "." + mineral.__name__
                    if not name in names:
                        names.add(name)
                        try:
                            x=mineral()
                            phases.append(x)
                        except:
                            print ""

        
    for a in names:
        print a


    print ""


    
    params = ['ref_V','ref_K','K_prime','ref_mu','mu_prime','molar_mass','n','ref_Debye','ref_grueneisen','q0','eta_0s']
    
    
    def create_list(mineral):
        row = [ p.__class__.__name__ ]
        for param in params:
            
            row.append(str(p.params[param] if param in p.params else ""))
        return row
    
    
    
    table = [ ['Name'] + params ]
    
    for p in phases:
        print p.__class__.__name__
        p.set_method('bm3')
        p.set_state(1e9,300)
        row = create_list(p)
        table.append(row)
            
    
    tools.pretty_print_table(table, False)
    
