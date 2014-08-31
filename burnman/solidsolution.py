# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

from burnman.mineral import Mineral
import warnings

import re
import numpy as np
from fractions import Fraction

R = 8.3145 # J/K/mol
kd = lambda x,y : 1 if x==y else 0

class SolidSolution(Mineral):
    """
    This is the base class for all solid solutions. 
    States of the solid solution can only be queried after setting 
    the pressure and temperature using set_state(). 
    The method for computing properties of
    the solution is set using set_method(), which should be done
    once after creating the material.

    This class is available as ``burnman.SolidSolution``.

    If deriving from this class, set the properties in self.params
    to the desired values. For more complicated materials you
    can overwrite set_state(), change the params and then call
    set_state() from this class.

    All the solid solution parameters are expected to be in SI units.  This
    means that the interaction parameters should be in J/mol, with the T 
    and P derivatives in J/K/mol and m^3/(mol molecule).
    """

    # init sets up matrices to speed up calculations for when P, T, X is defined.
    def __init__(self, base_material, interaction_parameter):
        # Initialise the solid soltion inputs
        self.base_material = base_material
        self.interaction_parameter = interaction_parameter

        # Number of endmembers in the solid solution
        self.n_endmembers=len(base_material)

        # Check that base_material and interaction parameter each have three parts
        assert(len(map(list, zip(*base_material))) == 3) 
        assert(len(interaction_parameter) == 3)
        
        # Check that the interaction parameters have the correct number of variables
        for i in range(3):
            assert(len(interaction_parameter[i]) == self.n_endmembers-1)
            for j in range(self.n_endmembers-1):
                assert(len(interaction_parameter[i][j]) == (self.n_endmembers-1)-j)

        # Split interaction parameter into H, S, V terms
        self.excess_enthalpy=interaction_parameter[0]
        self.excess_entropy=interaction_parameter[1]
        self.excess_volume=interaction_parameter[2]

        # Number of sites 
        self.n_sites=base_material[0][1].count('[')

        # Check the number of sites is the same for each endmember
        for i in range(self.n_endmembers):
            assert(base_material[i][1].count('[') == self.n_sites)

        # Check the chemical composition of each endmember is consistent with the dataset
        # NOT IMPLEMENTED YET

        # Number of unique site occupancies (e.g.. Mg on X etc.)
        self.sites=[[] for i in range(self.n_sites)]
        list_occupancies=[]
        list_multiplicity=np.empty(shape=(self.n_sites))
        self.n_occupancies=0
        for endmember in range(self.n_endmembers):
            list_occupancies.append([[0]*len(self.sites[site]) for site in range(self.n_sites)])
            s=re.split(r'\[', base_material[endmember][1])[1:]
            for site in range(self.n_sites):
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
            
                    if element_on_site not in self.sites[site]:
                        self.n_occupancies=self.n_occupancies+1
                        self.sites[site].append(element_on_site)
                        element_index=self.sites[site].index(element_on_site)
                        for parsed_mbr in range(len(list_occupancies)):
                            list_occupancies[parsed_mbr][site].append(0) 
                    else:
                        element_index=self.sites[site].index(element_on_site)
                    list_occupancies[endmember][site][element_index]=proportion_element_on_site

        # Site occupancies and multiplicities
        self.site_occupancies=np.empty(shape=(self.n_endmembers,self.n_occupancies))
        self.site_multiplicities=np.empty(shape=(self.n_occupancies))
        for endmember in range(self.n_endmembers):
            n_element=0
            for site in range(self.n_sites):
                for element in range(len(list_occupancies[endmember][site])):
                    self.site_occupancies[endmember][n_element]=list_occupancies[endmember][site][element]
                    self.site_multiplicities[n_element]=list_multiplicity[site]
                    n_element=n_element+1

        self.alpha=np.array([base_material[i][2] for i in range(self.n_endmembers)])

        self.Wh=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Ws=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Wv=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        for i in range(self.n_endmembers):
            for j in range(i+1, self.n_endmembers):
                self.Wh[i][j]=2.*self.excess_enthalpy[i][j-i-1]/(self.alpha[i]+self.alpha[j])
                self.Ws[i][j]=2.*self.excess_entropy[i][j-i-1]/(self.alpha[i]+self.alpha[j])
                self.Wv[i][j]=2.*self.excess_volume[i][j-i-1]/(self.alpha[i]+self.alpha[j])

        # END INITIALISATION


    def set_composition(self, molar_fraction):
        assert(len(self.base_material) == len(molar_fraction))
        assert(sum(molar_fraction) > 0.9999)
        assert(sum(molar_fraction) < 1.0001)

        self.molar_fraction=molar_fraction

        # Ideal activities
        self.occupancies=np.dot(self.molar_fraction, self.site_occupancies)
        self.ideal_activity=np.empty(shape=(self.n_endmembers))
        for endmember in range(self.n_endmembers):
            self.ideal_activity[endmember]=1.0
            normalisation_constant=1.0
            for element in range(self.n_occupancies):
                if self.site_occupancies[endmember][element] != 0:
                    self.ideal_activity[endmember]=self.ideal_activity[endmember]*pow(self.occupancies[element],self.site_multiplicities[element])
                    normalisation_constant=normalisation_constant/pow(self.site_occupancies[endmember][element],self.site_multiplicities[element])
            self.ideal_activity[endmember]=normalisation_constant*self.ideal_activity[endmember]

        # Nonideal contributions
        phi=np.array([self.alpha[i]*molar_fraction[i] for i in range(self.n_endmembers)])
        phi=np.divide(phi, np.sum(phi))
        self.H_excess=np.dot(self.alpha.T,molar_fraction)*np.dot(phi.T,np.dot(self.Wh,phi))
        self.S_excess=np.dot(self.alpha.T,molar_fraction)*np.dot(phi.T,np.dot(self.Ws,phi))
        self.V_excess=np.dot(self.alpha.T,molar_fraction)*np.dot(phi.T,np.dot(self.Wv,phi))

    def set_state(self, pressure, temperature, molar_fraction):
        # Set the state of all the endmembers
        for i in range(self.n_endmembers):
            self.base_material[i][0].set_method(self.base_material[i][0].params['equation_of_state'])
            self.base_material[i][0].set_state(pressure, temperature)

        # Find excess properties
        self.set_composition(molar_fraction)

        # Ideal contribution (configurational entropy)
        # Don't include endmembers with negligable molar fractions as the log term blows up
        tol=1e-10
        self.gibbs_excess_ideal=0
        for i in range(self.n_endmembers):
            if molar_fraction[i] > tol:
                self.gibbs_excess_ideal=self.gibbs_excess_ideal + molar_fraction[i]*R*temperature*np.log(self.ideal_activity[i])

        # Non-ideal contribution
        self.gibbs_excess_nonideal=self.H_excess - temperature*self.S_excess + pressure*self.V_excess

        # Total excess Gibbs
        self.gibbs_excess=self.gibbs_excess_ideal + self.gibbs_excess_nonideal

        self.V= sum([ self.base_material[i][0].V * self.molar_fraction[i] for i in range(self.n_endmembers) ]) + self.V_excess
        self.gibbs= sum([ self.base_material[i][0].gibbs * self.molar_fraction[i] for i in range(self.n_endmembers) ]) + self.gibbs_excess

        
        '''
        for prop in self.base_materials[0].params:
           try:
               self.params[prop] = sum([ self.base_materials[i].params[prop] * self.molar_fraction[i] for i in itrange ])
           except TypeError:
               #if there is a type error, it is probably a string.  Just go with the value of the first base_material.
               self.params[prop] = self.base_materials[0].params[prop]
        Mineral.set_state(self, pressure, temperature)
        '''


