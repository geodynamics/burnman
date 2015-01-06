# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import warnings

from burnman import Mineral
from solutionmodel import *


class SolidSolution(Mineral):
    """
    This is the base class for all solid solutions. 
    Site occupancies, endmember activities and the constant
    and pressure and temperature dependencies of the excess 
    properties can be queried after using set_composition()
    States of the solid solution can only be queried after setting 
    the pressure, temperature and composition using set_state().

    This class is available as :class:`burnman.SolidSolution`.
    It uses an instance of :class:`burnman.SolutionModel` to 
    calculate interaction terms between endmembers.

    All the solid solution parameters are expected to be in SI units.  This
    means that the interaction parameters should be in J/mol, with the T 
    and P derivatives in J/K/mol and m^3/mol.
    """

    def __init__(self, molar_fractions=None):
        """
        Set up matrices to speed up calculations for when P, T, X is defined.

        Parameters
        ----------
        endmembers: list of :class:`burnman.Mineral`
            List of endmembers in this solid solution.
        solution_model: :class:`burnman.SolutionModel`
            SolutionModel to use.
        """

        if hasattr(self, 'endmembers') == False:
            raise Exception("'endmembers' attribute missing from solid solution")

        # Set default solution model type 
        if hasattr(self, 'type'):
            if self.type == 'ideal':
                self.solution_model=IdealSolution(self.endmembers)
            else:
                if hasattr(self, 'enthalpy_interaction') == False:
                    self.enthalpy_interaction = None
                if hasattr(self, 'volume_interaction') == False:
                    self.volume_interaction = None
                if hasattr(self, 'entropy_interaction') == False:
                    self.entropy_interaction = None

                if self.type == 'symmetric':
                    self.solution_model=SymmetricRegularSolution(self.endmembers, self.enthalpy_interaction, self.volume_interaction, self.entropy_interaction)
                elif self.type == 'asymmetric':
                    try:
                        self.solution_model=AsymmetricRegularSolution(self.endmembers, self.alphas, self.enthalpy_interaction, self.volume_interaction, self.entropy_interaction)
                    except: 
                        raise Exception("'alphas' attribute missing from solid solution")
                elif self.type == 'subregular':
                    self.solution_model=SubregularSolution(self.endmembers, self.enthalpy_interaction, self.volume_interaction, self.entropy_interaction)
                else:
                    raise Exception("Solution model type "+self.params['type']+"not recognised.")
        else:
            warnings.warn("Warning, you have not set a solution model 'type' attribute for this solid solution.", stacklevel=2)
            self.solution_model=SolutionModel()

        # Number of endmembers in the solid solution
        self.n_endmembers = len(self.endmembers)

        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_method(self.endmembers[i][0].params['equation_of_state'])

        if molar_fractions is not None:
            self.molar_fractions = molar_fractions
            self.set_composition(molar_fractions)

    def get_endmembers(self):
        return self.endmembers

    def set_composition(self, molar_fractions ):
        """
        Set the composition for this solid solution.

        Parameters
        ----------
        molar_fractions: list of float
            molar abundance for each endmember, needs to sum to one.
        """
        assert(len(self.endmembers) == len(molar_fractions))
        assert(sum(molar_fractions) > 0.9999)
        assert(sum(molar_fractions) < 1.0001)
        self.molar_fractions = molar_fractions 

    def set_method(self, method):
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_method(method)
        self.method = self.endmembers[0][0].method

    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        """
        molar_mass = sum([ self.endmembers[i][0].molar_mass()*self.molar_fractions[i] for i in range(self.n_endmembers) ])
        return molar_mass

    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        """
        molar_mass = sum([ self.base_material[i][0].molar_mass()*self.molar_fraction[i] for i in range(self.n_endmembers) ])
        return molar_mass

    def set_state(self, pressure, temperature):
        self.pressure=pressure
        self.temperature=temperature
        # Set the state of all the endmembers
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_state(pressure, temperature)

        self.excess_partial_gibbs = self.solution_model.excess_partial_gibbs_free_energies( pressure, temperature, self.molar_fractions)
        self.excess_gibbs = self.solution_model.excess_gibbs_free_energy( pressure, temperature, self.molar_fractions)
        self.partial_gibbs = np.array([self.endmembers[i][0].gibbs for i in range(self.n_endmembers)]) + self.excess_partial_gibbs
        self.gibbs= sum([ self.endmembers[i][0].gibbs * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_gibbs

        self.excess_enthalpy = self.solution_model.excess_enthalpy( pressure, temperature, self.molar_fractions)
        self.excess_entropy = self.solution_model.excess_entropy( pressure, temperature, self.molar_fractions)
        self.excess_volume = self.solution_model.excess_volume( pressure, temperature, self.molar_fractions)

        self.H = sum([ self.endmembers[i][0].H * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_enthalpy
        self.S = sum([ self.endmembers[i][0].S * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_entropy
        self.V = sum([ self.endmembers[i][0].V * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_volume
        self.C_p = sum([ self.endmembers[i][0].C_p * self.molar_fractions[i] for i in range(self.n_endmembers) ])
        self.alpha = (1./self.V) * sum([ self.endmembers[i][0].alpha * self.endmembers[i][0].V * self.molar_fractions[i] for i in range(self.n_endmembers) ])
        self.K_T = self.V * 1./(sum([ self.endmembers[i][0].V / (self.endmembers[i][0].K_T)  * self.molar_fractions[i] for i in range(self.n_endmembers) ]))
 
        G_list = [ self.endmembers[i][0].G for i in range(self.n_endmembers) ]
        if 0.0 in G_list:
            self.G = 0.0
        else:
            self.G = self.V * 1./(sum([ self.endmembers[i][0].V / (self.endmembers[i][0].G)  * self.molar_fractions[i] for i in range(self.n_endmembers) ]))

        # Derived properties
        self.C_v = self.C_p - self.V*temperature*self.alpha*self.alpha*self.K_T

        # C_v and C_p -> 0 as T -> 0
        if temperature<1e-10:
            self.K_S = self.K_T
            self.gr = float('nan')
        else:
            self.K_S = self.K_T*self.C_p/self.C_v
            self.gr = self.alpha*self.K_T*self.V/self.C_v     

    def calcgibbs(self, pressure, temperature, molar_fractions): 
        return sum([ self.endmembers[i][0].calcgibbs(pressure, temperature) * molar_fractions[i] for i in range(self.n_endmembers) ]) + self.solution_model.excess_gibbs_free_energy( pressure, temperature, molar_fractions)

    def calcpartialgibbsexcesses(self, pressure, temperature, molar_fractions):
        return self.solution_model.excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions)
