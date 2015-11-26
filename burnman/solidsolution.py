# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


from __future__ import absolute_import

import numpy as np
import warnings

from .mineral import Mineral
from .solutionmodel import *


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

        # Endmember names
        self.endmember_names = []
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_method(self.endmembers[i][0].params['equation_of_state'])
            try:
                self.endmember_names.append(self.endmembers[i][0].params['name'])
            except AttributeError:
                self.endmember_names.append('')

        # Equation of state
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_method(self.endmembers[i][0].params['equation_of_state'])

        # Endmember compositions
        self.endmember_compositions=[self.endmembers[i][0].composition() for i in range(self.n_endmembers)]
        
        # Molar fractions
        if molar_fractions is not None:
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

        # Mineral composition
    def composition(self):
        bulk_composition=dict()
        for i, mbr_composition in enumerate(self.endmember_compositions):
            for element in mbr_composition:
                if element not in bulk_composition:
                    bulk_composition[element] = self.molar_fractions[i]*mbr_composition[element]
                else:
                    bulk_composition[element] += self.molar_fractions[i]*mbr_composition[element]
        return bulk_composition
                    
    def set_method(self, method):
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_method(method)
        self.method = self.endmembers[0][0].method

    '''
    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        """
        molar_mass = sum([ self.endmembers[i][0].molar_mass()*self.molar_fractions[i] for i in range(self.n_endmembers) ])
        return molar_mass
    '''

    def set_state(self, pressure, temperature):
        self.pressure=pressure
        self.temperature=temperature
        # Set the state of all the endmembers
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_state(pressure, temperature)
        
        '''
        self.excess_partial_gibbs = self.solution_model.excess_partial_gibbs_free_energies( pressure, temperature, self.molar_fractions)
        self.excess_gibbs = self.solution_model.excess_gibbs_free_energy( pressure, temperature, self.molar_fractions)
        self.partial_gibbs = np.array([self.endmembers[i][0].gibbs for i in range(self.n_endmembers)]) + self.excess_partial_gibbs
        #self.gibbs= sum([ self.endmembers[i][0].gibbs * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_gibbs
        
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
        '''

    def calcgibbs(self, pressure, temperature, molar_fractions): 
        return sum([ self.endmembers[i][0].calcgibbs(pressure, temperature) * molar_fractions[i] for i in range(self.n_endmembers) ]) + self.solution_model.excess_gibbs_free_energy( pressure, temperature, molar_fractions)

    def calcpartialgibbsexcesses(self, pressure, temperature, molar_fractions):
        return self.solution_model.excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions)


    @property
    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        """
        self._molar_mass = sum([ self.endmembers[i][0].molar_mass*self.molar_fractions[i] for i in range(self.n_endmembers) ])
        return self._molar_mass
    
    @property
    def density(self):
        """
        Returns density of the mineral [kg/m^3]
        """
        self._density = self.molar_mass/self.molar_volume
        return self._density
    
    @property
    def molar_volume(self):
        """
        Returns molar volume of the mineral [m^3/mol]
        """
        self.excess_volume = self.solution_model.excess_volume(self.pressure, self.temperature, self.molar_fractions)
        self._molar_volume=sum([ self.endmembers[i][0].molar_volume * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_volume
        return self._molar_volume
    
    @property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the mineral [unitless]
        """
        if self.temperature<1e-10:
            self._grueneisen_parameter = float('nan')
        else:
            self._grueneisen_parameter = self.thermal_expansivity*self.isothermal_bulk_modulus*self.molar_volume/self.heat_capacity_v
        return self._grueneisen_parameter
    
    @property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the mineral [Pa]
        """
        self._isothermal_bulk_modulus = self.V * 1./(sum([ self.endmembers[i][0].V / (self.endmembers[i][0].K_T)  * self.molar_fractions[i] for i in range(self.n_endmembers) ]))
        return self._isothermal_bulk_modulus
    
    @property
    def compressibility(self):
        """
        Returns compressibility of the mineral (or inverse isothermal bulk modulus) [1/Pa]
        """
        self._compressibility = 1./self.isothermal_bulk_modulus
        return self._compressibily
    
    @property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        """
        if self.temperature<1e-10:
            self._adiabatic_bulk_modulus = self.isothermal_bulk_modulus
        else:
            self._adiabatic_bulk_modulus = self.isothermal_bulk_modulus*self.heat_capacity_p/self.heat_capacity_v
        return self._adiabatic_bulk_modulus
    
    @property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        """
        G_list = [ self.endmembers[i][0].G for i in range(self.n_endmembers) ]
        if 0.0 in G_list:
            self._shear_modulus = 0.0
        else:
            self._shear_modulus = self.V * 1./(sum([ self.endmembers[i][0].V / (self.endmembers[i][0].G)  * self.molar_fractions[i] for i in range(self.n_endmembers) ]))
        return self._shear_modulus
    
    @property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha) of the mineral [1/K]
        """
        self._thermal_expansivity = (1./self.V) * sum([ self.endmembers[i][0].alpha * self.endmembers[i][0].V * self.molar_fractions[i] for i in range(self.n_endmembers) ])
        return self._thermal_expansivity
    
    @property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the mineral [J/K/mol]
        """
        self._heat_capacity_v = self.heat_capacity_p - self.molar_volume*self.temperature*self.thermal_expansivity*self.thermal_expansivity*self.isothermal_bulk_modulus
        return self._heat_capacity_v
    
    @property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the mineral [J/K/mol]
        """
        self._heat_capacity_p = sum([ self.endmembers[i][0].heat_capacity_p * self.molar_fractions[i] for i in range(self.n_endmembers) ])
        return self._heat_capacity_p
    
    @property
    def v_s(self):
        """
        Returns shear wave speed of the mineral [m/s]
        """
        self._v_s = np.sqrt(self.shear_modulus / self.density)
        return self._v_s
    
    @property
    def v_p(self):
        """
        Returns P wave speed of the mineral [m/s]
        """
        self._v_p = np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. * \
                             self.shear_modulus) / self.density)
        return self._v_p
    
    @property
    def v_phi(self):
        """
        Returns bulk sound speed of the mineral [m/s]
        """
        self._v_phi = np.sqrt(self.adiabatic_bulk_modulus / self.density)
        return self._v_phi
    
    
    @property
    def excess_partial_gibbs(self):
        """
        Returns excess partial gibbs free energy [J]
        Property specific to solid solutions.
        """
        self._excess_partial_gibbs = self.solution_model.excess_partial_gibbs_free_energies( self.pressure, self.temperature, self.molar_fractions)
        return self._excess_partial_gibbs

    @property
    def partial_gibbs(self):
        """
        Returns excess partial gibbs free energy [J]
        Property specific to solid solutions.
        """
        self._partial_gibbs = np.array([self.endmembers[i][0].gibbs for i in range(self.n_endmembers)]) + self.excess_partial_gibbs
        return self._partial_gibbs

 
    @property
    def excess_gibbs(self):
        """
        Returns excess gibbs free energy [J]
        Property specific to solid solutions.
        """
        self._excess_gibbs = self.solution_model.excess_gibbs_free_energy( self.pressure, self.temperature, self.molar_fractions)
        return self._excess_gibbs
    
    @property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the mineral [J]
        """
        self._gibbs= sum([ self.endmembers[i][0].gibbs * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_gibbs
        return self._gibbs
    
    @property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the mineral [J]
        """
        return None

    @property
    def excess_enthalpy(self):
        """
        Returns excess enthalpy [J]
        Property specific to solid solutions.
        """
        self._excess_enthalpy = self.solution_model.excess_enthalpy( self.pressure, self.temperature, self.molar_fractions)
        return self._excess_enthalpy
    
    @property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        """
        self._molar_enthalpy = sum([ self.endmembers[i][0].H * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_enthalpy
        return self._molar_enthalpy

    @property
    def excess_entropy(self):
        """
        Returns excess entropy [J]
        Property specific to solid solutions.
        """
        self._excess_entropy = self.solution_model.excess_entropy( self.pressure, self.temperature, self.molar_fractions)
        return self._excess_entropy
    
    @property
    def molar_entropy(self):
        """
        Returns enthalpy of the mineral [J]
        """
        self._molar_entropy = sum([ self.endmembers[i][0].S * self.molar_fractions[i] for i in range(self.n_endmembers) ]) + self.excess_entropy
        return self._molar_entropy

    gibbs = molar_gibbs
    V = molar_volume
    H = molar_enthalpy
    S  = molar_entropy
    C_p = heat_capacity_p
    C_v = heat_capacity_v
    alpha = thermal_expansivity
    K_T = isothermal_bulk_modulus
    K_S = adiabatic_bulk_modulus
    gr = grueneisen_parameter
    G = shear_modulus

