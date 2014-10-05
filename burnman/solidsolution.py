# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
from burnman.mineral import Mineral
from burnman.processchemistry import ProcessSolidSolutionChemistry
from burnman.solutionmodel import SolutionModel
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
    def __init__(self, base_material, solution_model = SolutionModel() ):
        # Initialise the solid solution inputs
        self.base_material = base_material
        self.solution_model = solution_model

        # Number of endmembers in the solid solution
        self.n_endmembers=len(base_material)

        for i in range(self.n_endmembers):
            self.base_material[i][0].set_method(self.base_material[i][0].params['equation_of_state'])

    def set_composition( self, molar_fraction ):
        assert(len(self.base_material) == len(molar_fraction))
        assert(sum(molar_fraction) > 0.9999)
        assert(sum(molar_fraction) < 1.0001)
        self.molar_fraction = molar_fraction 

    def set_state(self, pressure, temperature):
        self.pressure=pressure
        self.temperature=temperature
        # Set the state of all the endmembers
        for i in range(self.n_endmembers):
            self.base_material[i][0].set_state(pressure, temperature)

        self.excess_gibbs = self.solution_model.excess_gibbs_free_energy( pressure, temperature, self.molar_fraction)
        self.excess_volume = self.solution_model.excess_volume( pressure, temperature, self.molar_fraction)

        self.V= sum([ self.base_material[i][0].V * self.molar_fraction[i] for i in range(self.n_endmembers) ]) + self.excess_volume
        self.gibbs= sum([ self.base_material[i][0].gibbs * self.molar_fraction[i] for i in range(self.n_endmembers) ]) + self.excess_gibbs

        
        '''
        for prop in self.base_materials[0].params:
           try:
               self.params[prop] = sum([ self.base_materials[i].params[prop] * self.molar_fraction[i] for i in itrange ])
           except TypeError:
               #if there is a type error, it is probably a string.  Just go with the value of the first base_material.
               self.params[prop] = self.base_materials[0].params[prop]
        Mineral.set_state(self, pressure, temperature)
        '''


    def calcgibbs(self, pressure, temperature, molar_fractions): 
        return sum([ self.base_material[i][0].calcgibbs(pressure, temperature) * molar_fractions[i] for i in range(self.n_endmembers) ]) + self.solution_model.excess_gibbs_free_energy( pressure, temperature, molar_fractions)

    def calcpartialgibbsexcesses(self, pressure, temperature, molar_fractions):
        Hint, Sint, Vint = self.solution_model.non_ideal_interactions(molar_fractions)
        partialgibbsexcesses=np.empty(len(molar_fractions))
        partialgibbsexcesses = np.array([0.+R*temperature*self.solution_model.ln_ideal_activities(molar_fractions)[i] + Hint[i] - temperature*Sint[i] + pressure*Vint[i] for i in range(self.n_endmembers)])
        return partialgibbsexcesses
