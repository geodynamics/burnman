# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

""" Generates a text table with mineral properties. Run 'python table.py latex' to write a tex version of the table to mytable.tex """


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
                if inspect.isclass(mineral) and mineral!=burnman.mineral and issubclass(mineral, burnman.mineral):
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


    
    params = ['V_0','K_0','Kprime_0','G_0','Gprime_0','molar_mass','n','Debye_0','grueneisen_0','q_0','eta_s_0']
    
    
    def create_list(mineral):
        row = [ p.to_string() ]
        for param in params:
            
            row.append(str(p.params[param] if param in p.params else ""))
        return row
    
    
    
    table = [ ['Name'] + params ]
    tablel = []

    for p in phases:
        print p.to_string()
        p.set_method('bm3')
        p.set_state(1e9,300)
        row = create_list(p)
        table.append(row)
        tablel.append(row)    
 
    if (len(sys.argv)==1):
        tools.pretty_print_table(table, False)
    elif sys.argv[1]=="tab":
        tools.pretty_print_table(table, True)
    elif sys.argv[1]=="latex":
        import table_latex as Table

        for r in tablel:
            r[0]=r[0].replace("'","").replace("_","\_")

        fout = open('mytable.tex','w')
        t = Table.Table(12, justs='lrrrrrrrrrrr', caption='Mineral parameters used', label="tab:label", rotate='True', fontsize='footnotesize')
        t.add_header_row(['Name','V\_0','K\_0','K\_prime','G\_0','G\_prime','molar\_mass','n','Debye\_0','grueneisen\_0','q\_0','eta\_0s'] )
        t.add_data(map(list,zip(*tablel)), sigfigs=2)
        t.print_table(fout)
        fout.close() 
        print "see mytable.tex for the latex version of the table" 
