# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# This module provides the functions required to process the standard burnman formula compositions
# ProcessChemistry returns the number of atoms and molar mass of a compound given its unit formula as an argument.
# ProcessSolidSolutionChemistry returns information required to calculate solid solution properties from a set of endmember formulae

import re
import numpy as np
from fractions import Fraction

def read_masses(filename): 
    lookup=dict()
    for line in open(filename, "r"):
        data="%".join(line.split("%")[:1]).split()
        if data != []:
            lookup[data[0]]=float(data[1])
    return lookup

def dictionarize_formula(formula):
    f=dict()
    for element in re.findall('[A-Z][^A-Z]*', formula):
        list=filter(None, re.split(r'(\d+)', element))
        # Look up number of atoms of element
        if len(list) == 1:
            nel=1.
        else: 
            nel=float(list[1])
        f[list[0]]=nel
    return f

def formula_mass(formula, atomic_masses):
    mass=0.
    for element in formula:
        mass=mass+formula[element]*atomic_masses[element]
    return mass

def dictionarize_site_formula(formula):
    solution_formulae=dict()
    sites=[[] for i in range(len(s))]
    list_occupancies=[]
    list_multiplicity=np.empty(shape=(len(s)))
    n_occupancies=0

    s=re.split(r'\[', formula)[1:]
    
    for site in range(len(s)):
        site_occupancy=re.split(r'\]', s[site])[0]
        mult=re.split('[A-Z][^A-Z]*',re.split(r'\]', s[site])[1])[0]
        not_in_site=filter(None, re.split(r'\]', s[site]))[1]
        not_in_site=not_in_site.replace(mult, '', 1)
        if mult == '':
            list_multiplicity[site]=1.0
        else:
            list_multiplicity[site]=mult

        # Loop over elements on a site
        elements=re.findall('[A-Z][^A-Z]*',site_occupancy)
        for i in range(len(elements)):
            element_on_site=re.split('[0-9][^A-Z]*',elements[i])[0]
            proportion_element_on_site=re.findall('[0-9][^A-Z]*',elements[i])
            if len(proportion_element_on_site) == 0:
                proportion_element_on_site=Fraction(1.0)
            else:
                proportion_element_on_site=Fraction(proportion_element_on_site[0])
            n_element=float(mult)*proportion_element_on_site
            f[element_on_site]=f.get(element_on_site, 0.0) + n_element

        # Loop over elements not on a site
        for enamenumber in re.findall('[A-Z][^A-Z]*', not_in_site):
            element=filter(None, re.split(r'(\d+)', enamenumber))
            f[element[0]]=f.get(element[0], 0.0) + float(element[1])

    return f


def ProcessSolidSolutionChemistry(formula):
    n_sites=formula[0].count('[')
    n_endmembers=len(formula)

    # Check the number of sites is the same for all endmembers
    for i in range(n_endmembers):
        assert(formula[i].count('[') == n_sites)

    solution_formulae=[]
    sites=[[] for i in range(n_sites)]
    list_occupancies=[]
    list_multiplicity=np.empty(shape=(n_sites))
    n_occupancies=0
        
    s=re.split(r'\[', formula[i])[1:]

        # Number of unique site occupancies (e.g.. Mg on X etc.)
    for endmember in range(n_endmembers):
        solution_formula=dict()
        list_occupancies.append([[0]*len(sites[site]) for site in range(n_sites)])
        s=re.split(r'\[', formula[endmember])[1:]

        for site in range(n_sites):
            site_occupancy=re.split(r'\]', s[site])[0]
            mult=re.split('[A-Z][^A-Z]*',re.split(r'\]', s[site])[1])[0]
            not_in_site=filter(None, re.split(r'\]', s[site]))[1]
            not_in_site=not_in_site.replace(mult, '', 1)
            if mult == '':
                list_multiplicity[site]=1.0
            else:
                list_multiplicity[site]=mult

            # Loop over elements on a site
            elements=re.findall('[A-Z][^A-Z]*',site_occupancy)
            for i in range(len(elements)):
                element_on_site=re.split('[0-9][^A-Z]*',elements[i])[0]
                proportion_element_on_site=re.findall('[0-9][^A-Z]*',elements[i])
                if len(proportion_element_on_site) == 0:
                    proportion_element_on_site=Fraction(1.0)
                else:
                    proportion_element_on_site=Fraction(proportion_element_on_site[0])
                solution_formula[element_on_site]=solution_formula.get(element_on_site, 0.0) + float(mult)*proportion_element_on_site
            
                if element_on_site not in sites[site]:
                    n_occupancies=n_occupancies+1
                    sites[site].append(element_on_site)
                    element_index=sites[site].index(element_on_site)
                    for parsed_mbr in range(len(list_occupancies)):
                        list_occupancies[parsed_mbr][site].append(0) 
                else:
                    element_index=sites[site].index(element_on_site)
                list_occupancies[endmember][site][element_index]=proportion_element_on_site

            # Loop over elements not on a site
            for enamenumber in re.findall('[A-Z][^A-Z]*', not_in_site):
                element=filter(None, re.split(r'(\d+)', enamenumber))
                solution_formula[element[0]]=solution_formula.get(element[0], 0.0) + float(element[1])

        solution_formulae.append(solution_formula)

    # Site occupancies and multiplicities
    endmember_occupancies=np.empty(shape=(n_endmembers,n_occupancies))
    site_multiplicities=np.empty(shape=(n_occupancies))
    for endmember in range(n_endmembers):
        n_element=0
        for site in range(n_sites):
            for element in range(len(list_occupancies[endmember][site])):
                endmember_occupancies[endmember][n_element]=list_occupancies[endmember][site][element]
                site_multiplicities[n_element]=list_multiplicity[site]
                n_element=n_element+1
    

    return solution_formulae, n_sites, sites, n_occupancies, endmember_occupancies, site_multiplicities

