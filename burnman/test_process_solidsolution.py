# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

# This is a function which returns the number of atoms and molar mass of a compound given its unit formula as an argument.

import re

        # Endmembers
base_material = [['pyrope()',    '[Mg]3[Al]2Si3O12',         1.0], \
                 ['almandine()', '[Fe]3[Al]2Si3O12',         1.0], \
                 ['grossular()', '[Ca]3[Al]2Si3O12',         2.7], \
                 ['majorite()',  '[Mg]3[Mg1/2Si1/2]2Si3O12', 1.0]]


nsites=base_material[0][1].count('[')
print 'Number of sites:', nsites

for k in range(len(base_material)):
    print base_material[k][0]
    sites=re.split(r'\[', base_material[k][1])[1:]
    for i in range(nsites):
        site_occupancy=re.split(r'\]', sites[i])[0]
        site_multiplicity=re.split('[A-Z][^A-Z]*',re.split(r'\]', sites[i])[1])[0]
        print 'Multiplicity on site', i+1, ':', site_multiplicity
        elements=re.findall('[A-Z][^A-Z]*',site_occupancy)
        for j in range(len(elements)):
            element_on_site=re.split('[0-9][^A-Z]*',elements[j])[0]
            proportion_element_on_site=re.findall('[0-9][^A-Z]*',elements[j])
            if len(proportion_element_on_site) == 0:
                proportion_element_on_site=1.0
            else:
                proportion_element_on_site=proportion_element_on_site[0]
            print element_on_site, proportion_element_on_site
    print ''

print ''


sites=[[] for i in range(nsites)]
for i in range(len(base_material)):
    print base_material[i][1]
    for j in range(nsites):
        sites[j].append(j)
print sites    

#2D vector of endmembers and alphas
#2D vector of sites and multiplicities [[A,3],[B,2]]
#2D vector of site occupancies e.g. [[x(Mg,A),x(Fe,A)],[x(Al,B),x(Mg,B),x(Si,B)]]
#3D vector of site occupancies [[[1.,0.],[1.,0.,0.]],[[0.,1.],[1.,0.,0.]]]
#2D vector of interaction parameters and temperature and pressure derivatives


#Operations required:
#- matrix operation to return site occupancies, given proportions of endmembers
#- operation to obtain ideal activities
#- operation to obtain non ideal activities
