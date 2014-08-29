# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# This is a test script which is the precursor to implementing a solid solution BaseClass.

import re
import numpy as np
from fractions import Fraction

# Endmembers
base_material = [['pyrope()',    '[Mg]3[Al]2Si3O12',         1.0], \
                 ['almandine()', '[Fe]3[Al]2Si3O12',         1.0], \
                 ['grossular()', '[Ca]3[Al]2Si3O12',         2.7], \
                 ['majorite()',  '[Mg]3[Mg1/2Si1/2]2Si3O12', 1.0]]

n_endmembers=len(base_material)

# Interaction parameters
excess_enthalpy=[[2.5e3, 29.1e3, 15e3],[10e3,18e3],[48e3]]
excess_entropy=[[0., 0., 0.],[0., 0.],[0.]]
excess_volume=[[0., 0., 0.],[0., 0.],[0.]]

interaction_parameter=[excess_enthalpy,excess_entropy,excess_volume]



endmember_proportions = np.array([ 0.5, 0.2, 0.1, 0.2 ])
#endmember_proportions = np.array([ 0.0, 0.0, 0.0, 1.0 ])
#endmember_proportions = np.array([ 1.0, 0.0, 0.0, 0.0 ])

# "sites" is a 2D list of sites and the elements which reside on them 
# "site_occupancies" is a 3D list describing the elemental site occupancies of each endmember 
# "array_occupancies" is a 2D np.array of site_occupancies, concatenating the 2nd and 3rd dimension
nsites=base_material[0][1].count('[')
print 'Number of sites:', nsites
print ''

sites=[[] for i in range(nsites)]
site_occupancies=[]
site_multiplicity=np.empty(shape=(nsites))
n_elements=0
for endmember in range(n_endmembers):
    site_occupancies.append([[0]*len(sites[site]) for site in range(nsites)])
    s=re.split(r'\[', base_material[endmember][1])[1:]
    for site in range(nsites):
        site_occupancy=re.split(r'\]', s[site])[0]
        mult=re.split('[A-Z][^A-Z]*',re.split(r'\]', s[site])[1])[0]
        if mult == '':
            site_multiplicity[site]=1.0
        else:
            site_multiplicity[site]=mult
        elements=re.findall('[A-Z][^A-Z]*',site_occupancy)
        for i in range(len(elements)):
            element_on_site=re.split('[0-9][^A-Z]*',elements[i])[0]
            proportion_element_on_site=re.findall('[0-9][^A-Z]*',elements[i])
            if len(proportion_element_on_site) == 0:
                proportion_element_on_site=Fraction(1.0)
            else:
                proportion_element_on_site=Fraction(proportion_element_on_site[0])
            
            if element_on_site not in sites[site]:
                n_elements=n_elements+1
                sites[site].append(element_on_site)
                element_index=sites[site].index(element_on_site)
                for parsed_mbr in range(len(site_occupancies)):
                    site_occupancies[parsed_mbr][site].append(0) 
            else:
                element_index=sites[site].index(element_on_site)
            site_occupancies[endmember][site][element_index]=proportion_element_on_site



array_occupancies=np.empty(shape=(n_endmembers,n_elements))
array_multiplicities=np.empty(shape=(n_elements))
for endmember in range(n_endmembers):
    n_element=0
    for site in range(nsites):
        for element in range(len(site_occupancies[endmember][site])):
            array_occupancies[endmember][n_element]=site_occupancies[endmember][site][element]
            array_multiplicities[n_element]=site_multiplicity[site]
            n_element=n_element+1


# Matrix operation to return site occupancies, given proportions of endmembers
occupancies=np.dot(endmember_proportions, array_occupancies)
print 'Site occupancies'
print sites
print occupancies
print ''



# Ideal activities
ideal_activity=np.empty(shape=(n_endmembers))
for endmember in range(n_endmembers):
    ideal_activity[endmember]=1.0
    normalisation_constant=1.0
    for element in range(n_elements):
        if array_occupancies[endmember][element] != 0:
            #print base_material[endmember][0], element, occupancies[element], array_occupancies[endmember][element]
            ideal_activity[endmember]=ideal_activity[endmember]*pow(occupancies[element],array_multiplicities[element])
            normalisation_constant=normalisation_constant/pow(array_occupancies[endmember][element],array_multiplicities[element])
    ideal_activity[endmember]=normalisation_constant*ideal_activity[endmember]


print 'Ideal activities'
print ideal_activity

# Non ideal activities
# alpha^T*p*(phi^T*W*phi)

alpha=np.array([base_material[i][2] for i in range(n_endmembers)])
phi=np.array([alpha[i]*endmember_proportions[i] for i in range(n_endmembers)])
phi=np.divide(phi, np.sum(phi))

Wh=np.zeros(shape=(n_endmembers,n_endmembers))
Ws=np.zeros(shape=(n_endmembers,n_endmembers))
Wv=np.zeros(shape=(n_endmembers,n_endmembers))
for i in range(n_endmembers):
    for j in range(i+1, n_endmembers):
        Wh[i][j]=2*excess_enthalpy[i][j-i-1]/(alpha[i]+alpha[j])
        Ws[i][j]=2*excess_entropy[i][j-i-1]/(alpha[i]+alpha[j])
        Wv[i][j]=2*excess_volume[i][j-i-1]/(alpha[i]+alpha[j])

print ''
print 'Nonideal contribution to gibbs'
print np.dot(alpha.T,endmember_proportions)*np.dot(phi.T,np.dot(Wh,phi))
