# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np

from burnman import Mineral
from solutionmodel import SolutionModel
from solutionmodel import kd


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

    def __init__(self, endmembers, solution_model=SolutionModel()):
        """
        Set up matrices to speed up calculations for when P, T, X is defined.

        Parameters
        ----------
        endmembers: list of :class:`burnman.Mineral`
            List of endmembers in this solid solution.
        solution_model: :class:`burnman.SolutionModel`
            SolutionModel to use.
        """

        # Initialise the solid solution inputs
        self.base_material = endmembers
        self.solution_model = solution_model

        # Number of endmembers in the solid solution
        self.n_endmembers = len(endmembers)

        for i in range(self.n_endmembers):
            self.base_material[i][0].set_method(self.base_material[i][0].params['equation_of_state'])

    def get_endmembers(self):
        return self.base_material

    def set_composition(self, molar_fraction ):
        """
        Set the composition for this solid solution.

        Parameters
        ----------
        molar_fraction: list of float
            molar abundance for each endmember, needs to sum to one.
        """
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
        self.alpha = (1./self.V) * sum([ self.base_material[i][0].alpha * self.base_material[i][0].V * self.molar_fraction[i] for i in range(self.n_endmembers) ])
        self.K_T = self.V * 1./(sum([ self.base_material[i][0].V / (self.base_material[i][0].K_T)  * self.molar_fraction[i] for i in range(self.n_endmembers) ]))
        self.G = self.V * 1./(sum([ self.base_material[i][0].V / (self.base_material[i][0].G)  * self.molar_fraction[i] for i in range(self.n_endmembers) ]))

        # Derived properties
        self.C_v = self.C_p - self.V*temperature*self.alpha*self.alpha*self.K_T

        # C_v and C_p -> 0 as T -> 0
        if temperature<1e-10:
            self.K_S = self.K_T
            self.gr = float('nan')
        else:
            self.K_S = self.K_T*self.C_p/self.C_v
            self.gr = self.alpha*self.K_T*self.V/self.C_v     

    def calcgibbs(self, pressure, temperature, molar_fraction): 
        return sum([ self.base_material[i][0].calcgibbs(pressure, temperature) * molar_fraction[i] for i in range(self.n_endmembers) ]) + self.solution_model.excess_gibbs_free_energy( pressure, temperature, molar_fraction)

    def calcpartialgibbsexcesses(self, pressure, temperature, molar_fraction):
        return self.solution_model.excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fraction)
