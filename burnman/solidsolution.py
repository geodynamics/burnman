# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
from burnman.mineral import Mineral
from burnman.processchemistry import ProcessSolidSolutionChemistry
from burnman.solutionmodel import SolutionModel
import constants
import warnings

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

        self.excess_partial_gibbs = self.solution_model.excess_partial_gibbs_free_energies( pressure, temperature, self.molar_fraction)
        self.excess_gibbs = self.solution_model.excess_gibbs_free_energy( pressure, temperature, self.molar_fraction)
        self.partial_gibbs = np.array([self.base_material[i][0].gibbs for i in range(self.n_endmembers)]) + self.excess_partial_gibbs
        self.gibbs= sum([ self.base_material[i][0].gibbs * self.molar_fraction[i] for i in range(self.n_endmembers) ]) + self.excess_gibbs

        self.excess_enthalpy = self.solution_model.excess_enthalpy( pressure, temperature, self.molar_fraction)
        self.excess_entropy = self.solution_model.excess_entropy( pressure, temperature, self.molar_fraction)
        self.excess_volume = self.solution_model.excess_volume( pressure, temperature, self.molar_fraction)

        self.H = sum([ self.base_material[i][0].H * self.molar_fraction[i] for i in range(self.n_endmembers) ]) + self.excess_enthalpy
        self.S = sum([ self.base_material[i][0].S * self.molar_fraction[i] for i in range(self.n_endmembers) ]) + self.excess_entropy
        self.V = sum([ self.base_material[i][0].V * self.molar_fraction[i] for i in range(self.n_endmembers) ]) + self.excess_volume
        self.C_p = sum([ self.base_material[i][0].C_p * self.molar_fraction[i] for i in range(self.n_endmembers) ])

        # The following are currently just simple molar sums ...
        # ... when they are corrected, they will be moved up ...
        self.gr = sum([ self.base_material[i][0].gr * self.molar_fraction[i] for i in range(self.n_endmembers) ])
        self.K_T = sum([ self.base_material[i][0].K_T * self.molar_fraction[i] for i in range(self.n_endmembers) ])
        self.K_S = sum([ self.base_material[i][0].K_S * self.molar_fraction[i] for i in range(self.n_endmembers) ])
        self.C_v = sum([ self.base_material[i][0].C_v * self.molar_fraction[i] for i in range(self.n_endmembers) ])
        self.alpha = sum([ self.base_material[i][0].alpha * self.molar_fraction[i] for i in range(self.n_endmembers) ])


    def calcgibbs(self, pressure, temperature, molar_fractions): 
        return sum([ self.base_material[i][0].calcgibbs(pressure, temperature) * molar_fractions[i] for i in range(self.n_endmembers) ]) + self.solution_model.excess_gibbs_free_energy( pressure, temperature, molar_fractions)

    def calcpartialgibbsexcesses(self, pressure, temperature, molar_fractions):
        return self.solution_model.excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions)
