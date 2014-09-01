# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
from burnman.mineral import Mineral
from burnman.processchemistry import ProcessSolidSolutionChemistry
from burnman.processchemistry import CompositionEquality
import warnings

R = 8.3145 # J/K/mol
kd = lambda x,y : 1 if x==y else 0

class SolidSolution(Mineral):
    """
    This is the base class for all solid solutions. 
    Site occupancies, endmember activities and the constant
    and pressure and temperature dependencies of the excess 
    properties can be queried after using set_composition()
    States of the solid solution can only be queried after setting 
    the pressure, temperature and composition using set_state().

    This class is available as ``burnman.SolidSolution``.

    All the solid solution parameters are expected to be in SI units.  This
    means that the interaction parameters should be in J/mol, with the T 
    and P derivatives in J/K/mol and m^3/mol.
    """

    # init sets up matrices to speed up calculations for when P, T, X is defined.
    def __init__(self, base_material, interaction_parameter):
        # Initialise the solid solution inputs
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

        # Check the chemical composition of each endmember is consistent with the dataset
        for i in range(self.n_endmembers):
            endmember_formula=base_material[i][0].params['formula']
            solution_formula=base_material[i][1]
            if CompositionEquality(endmember_formula, solution_formula) != True:
                print 'Formula of endmember', base_material[i][0].params['name'], 'does not agree with formula in the', self.name, 'SolidSolution model'
                exit()

        # Process solid solution chemistry
        self.n_sites, self.sites, self.n_occupancies, self.endmember_occupancies, self.site_multiplicities = ProcessSolidSolutionChemistry([base_material[i][1] for i in range(self.n_endmembers)])

        # Create array of van Laar parameters
        self.alpha=np.array([base_material[i][2] for i in range(self.n_endmembers)])

        # Create 2D arrays of interaction parameters
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
        # Check that the composition declared is consistent with the solution model
        assert(len(self.base_material) == len(molar_fraction))
        assert(sum(molar_fraction) > 0.9999)
        assert(sum(molar_fraction) < 1.0001)

        # Make molar fraction an attribute of the solid solution
        self.molar_fraction=molar_fraction

        # Ideal activities
        self.site_occupancies=np.dot(self.molar_fraction, self.endmember_occupancies)
        self.ideal_activity=np.empty(shape=(self.n_endmembers))
        for endmember in range(self.n_endmembers):
            self.ideal_activity[endmember]=1.0
            normalisation_constant=1.0
            for element in range(self.n_occupancies):
                if self.endmember_occupancies[endmember][element] != 0:
                    self.ideal_activity[endmember]=self.ideal_activity[endmember]*pow(self.site_occupancies[element],self.site_multiplicities[element])
                    normalisation_constant=normalisation_constant/pow(self.endmember_occupancies[endmember][element],self.site_multiplicities[element])
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
