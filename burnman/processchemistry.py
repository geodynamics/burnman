# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# This is a function which returns the number of atoms and molar mass of a compound given its unit formula as an argument.

import re

def ProcessChemistry(formula):
    filename="data/input_masses/atomic_masses.dat"
    file = open(filename, "r")
    el_name = []
    el_mass = []
    
    i=0
    for line in file:
        data="%".join(line.split("%")[:1]).split()
        if data != []:
            el_name.append(data[0])
            el_mass.append(float(data[1]))
    file.close()

# Loop over elements
    n=0.
    molar_mass=0.
    for element in re.findall('[A-Z][^A-Z]*', formula):
        list=filter(None, re.split(r'(\d+)', element))
    # Look up number of atoms of element
        if len(list) == 1:
            nel=1.
        else: 
            nel=float(list[1])
    # Increment size of compound
        n=n+nel

    # Find atomic mass of element
        molar_mass=molar_mass+nel*el_mass[el_name.index(list[0])]

    return n, molar_mass
