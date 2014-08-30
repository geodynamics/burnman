# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# This is a test script which is the precursor to implementing a solid solution BaseClass.


'''
# Inputs
Solid solution model
Endmember proportions
P,T

# Things we can initialise without endmember proportions
n_sites
n_occupancies
n_endmembers
alpha
Wh, Ws, Wv
site_occupancies
site_multiplicities

# Endmember proportions needed
phi
ideal activities
occupancies

# P, T needed
Wtotal / nonideal contribution to gibbs
'''

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
excess_enthalpy=[[2.5e3, 30.1e3, 15e3],[10e3,18e3],[48e3]]
excess_entropy=[[0., 0., 0.],[0., 0.],[0.]]
excess_volume=[[0., 0.164e-5, 0.],[0., 0.],[0.]]
interaction_parameter=[excess_enthalpy,excess_entropy,excess_volume]


# INPUT PROPORTIONS
molar_fraction = np.array([ 0.5, 0.2, 0.1, 0.2 ])


# "sites" is a 2D list of sites and the elements which reside on them 
# "list_occupancies" is a 3D list describing the elemental site occupancies of each endmember 
# "site_occupancies" is a 2D np.array of list_occupancies, concatenating the 2nd and 3rd dimension
n_sites=base_material[0][1].count('[')
print 'Number of sites:', n_sites
print ''

sites=[[] for i in range(n_sites)]
list_occupancies=[]
list_multiplicity=np.empty(shape=(n_sites))
n_occupancies=0
for endmember in range(n_endmembers):
    list_occupancies.append([[0]*len(sites[site]) for site in range(n_sites)])
    s=re.split(r'\[', base_material[endmember][1])[1:]
    for site in range(n_sites):
        site_occupancy=re.split(r'\]', s[site])[0]
        mult=re.split('[A-Z][^A-Z]*',re.split(r'\]', s[site])[1])[0]
        if mult == '':
            list_multiplicity[site]=1.0
        else:
            list_multiplicity[site]=mult
        elements=re.findall('[A-Z][^A-Z]*',site_occupancy)
        for i in range(len(elements)):
            element_on_site=re.split('[0-9][^A-Z]*',elements[i])[0]
            proportion_element_on_site=re.findall('[0-9][^A-Z]*',elements[i])
            if len(proportion_element_on_site) == 0:
                proportion_element_on_site=Fraction(1.0)
            else:
                proportion_element_on_site=Fraction(proportion_element_on_site[0])
            
            if element_on_site not in sites[site]:
                n_occupancies=n_occupancies+1
                sites[site].append(element_on_site)
                element_index=sites[site].index(element_on_site)
                for parsed_mbr in range(len(list_occupancies)):
                    list_occupancies[parsed_mbr][site].append(0) 
            else:
                element_index=sites[site].index(element_on_site)
            list_occupancies[endmember][site][element_index]=proportion_element_on_site



site_occupancies=np.empty(shape=(n_endmembers,n_occupancies))
site_multiplicities=np.empty(shape=(n_occupancies))
for endmember in range(n_endmembers):
    n_element=0
    for site in range(n_sites):
        for element in range(len(list_occupancies[endmember][site])):
            site_occupancies[endmember][n_element]=list_occupancies[endmember][site][element]
            site_multiplicities[n_element]=list_multiplicity[site]
            n_element=n_element+1


# Matrix operation to return site occupancies, given proportions of endmembers
occupancies=np.dot(molar_fraction, site_occupancies)
print 'Site occupancies'
print sites
print occupancies
print ''



# Ideal activities
ideal_activity=np.empty(shape=(n_endmembers))
for endmember in range(n_endmembers):
    ideal_activity[endmember]=1.0
    normalisation_constant=1.0
    for element in range(n_occupancies):
        if site_occupancies[endmember][element] != 0:
            #print base_material[endmember][0], element, occupancies[element], site_occupancies[endmember][element]
            ideal_activity[endmember]=ideal_activity[endmember]*pow(occupancies[element],site_multiplicities[element])
            normalisation_constant=normalisation_constant/pow(site_occupancies[endmember][element],site_multiplicities[element])
    ideal_activity[endmember]=normalisation_constant*ideal_activity[endmember]


print 'Ideal activities'
print ideal_activity

# Non ideal activities
# alpha^T*p*(phi^T*W*phi)

alpha=np.array([base_material[i][2] for i in range(n_endmembers)])
phi=np.array([alpha[i]*molar_fraction[i] for i in range(n_endmembers)])
phi=np.divide(phi, np.sum(phi))

Wh=np.zeros(shape=(n_endmembers,n_endmembers))
Ws=np.zeros(shape=(n_endmembers,n_endmembers))
Wv=np.zeros(shape=(n_endmembers,n_endmembers))
for i in range(n_endmembers):
    for j in range(i+1, n_endmembers):
        Wh[i][j]=2.*excess_enthalpy[i][j-i-1]/(alpha[i]+alpha[j])
        Ws[i][j]=2.*excess_entropy[i][j-i-1]/(alpha[i]+alpha[j])
        Wv[i][j]=2.*excess_volume[i][j-i-1]/(alpha[i]+alpha[j])

print ''
print 'Excess enthalpy, entropy, volume (non-configurational)'
print np.dot(alpha.T,molar_fraction)*np.dot(phi.T,np.dot(Wh,phi)), 'J/mol'

print np.dot(alpha.T,molar_fraction)*np.dot(phi.T,np.dot(Ws,phi)), 'J/K/mol'

print np.dot(alpha.T,molar_fraction)*np.dot(phi.T,np.dot(Wv,phi)), 'm^3/mol'



# Plot excess volumes for the pyrope-grossular join
ppy= np.linspace(0, 1, 101)
vex= np.empty(shape=(101))
for p in range(len(ppy)):
    a=ppy[p]
    molar_fraction = np.array([ a, 0.0, 1.-a, 0.0 ])
    phi=np.array([alpha[i]*molar_fraction[i] for i in range(n_endmembers)])
    phi=np.divide(phi, np.sum(phi))

    vex[p]=np.dot(alpha.T,molar_fraction)*np.dot(phi.T,np.dot(Wv,phi))


import matplotlib.pyplot as plt
plt.plot(ppy,vex,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=0)
plt.xlim(min(ppy),max(ppy))
plt.xlabel("p(pyrope)")
plt.title("V excess (m^3/mol) for pyrope-grossular garnets")

plt.show()
